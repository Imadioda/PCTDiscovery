from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone
from urllib.request import urlopen
import hashlib
import json

import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "rag_data"

INPUT_FILE = DATA_DIR / "retained_genes.csv"
REFERENCE_DIR = DATA_DIR / "references"
HGNC_FILE = REFERENCE_DIR / "hgnc_complete_set.txt"

HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/"
    "hgnc/tsv/tsv/hgnc_complete_set.txt"
)


def load_reference():
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

    if not HGNC_FILE.exists():
        print("Downloading the HGNC reference...")

        with urlopen(HGNC_URL, timeout=120) as response:
            content = response.read()

        # Validate before storing the download.
        from io import BytesIO
        reference = pd.read_csv(
            BytesIO(content), sep="\t", dtype=str, keep_default_na=False
        )
        validate_reference(reference)

        temporary = HGNC_FILE.with_suffix(".tmp")
        temporary.write_bytes(content)
        temporary.replace(HGNC_FILE)

        provenance = {
            "source_url": HGNC_URL,
            "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
            "sha256": hashlib.sha256(content).hexdigest(),
        }
        HGNC_FILE.with_suffix(".metadata.json").write_text(
            json.dumps(provenance, indent=2), encoding="utf-8"
        )

    reference = pd.read_csv(
        HGNC_FILE, sep="\t", dtype=str, keep_default_na=False
    )
    validate_reference(reference)
    return reference.loc[reference["status"].eq("Approved")].copy()


def validate_reference(reference):
    required = {
        "hgnc_id", "symbol", "name", "status",
        "prev_symbol", "alias_symbol", "ensembl_gene_id",
    }
    missing = required - set(reference.columns)

    if missing or reference.empty:
        raise ValueError(
            f"Invalid HGNC reference. Missing columns: {sorted(missing)}"
        )

    if reference["hgnc_id"].duplicated().any():
        raise ValueError("Duplicate HGNC identifiers in the reference.")


def build_indexes(reference):
    current = defaultdict(set)
    historical = defaultdict(set)

    for row in reference.to_dict("records"):
        identifier = row["hgnc_id"]
        current[row["symbol"]].add(identifier)

        for field in ["prev_symbol", "alias_symbol"]:
            for label in row[field].split("|"):
                label = label.strip()
                if label:
                    historical[label].add(identifier)

    records = reference.set_index("hgnc_id").to_dict("index")
    return current, historical, records


def annotate_gene(gene, current, historical, records):
    label = gene.strip()

    # Exact official symbols take precedence.
    candidates = current.get(label, set())
    method = "approved_symbol"

    if not candidates:
        candidates = historical.get(label, set())
        method = "previous_or_alias_symbol"

    annotation = {
        "gene": gene,
        "mapping_status": "unresolved",
        "mapping_method": "",
        "hgnc_id": "",
        "approved_symbol": "",
        "gene_name": "",
        "ensembl_gene_id": "",
        "previous_symbols": "",
        "alias_symbols": "",
        "candidate_hgnc_ids": "|".join(sorted(candidates)),
    }

    if not candidates:
        return annotation

    if len(candidates) > 1:
        annotation["mapping_status"] = "ambiguous"
        return annotation

    identifier = next(iter(candidates))
    record = records[identifier]

    annotation.update({
        "mapping_status": "mapped",
        "mapping_method": method,
        "hgnc_id": identifier,
        "approved_symbol": record["symbol"],
        "gene_name": record["name"],
        "ensembl_gene_id": record["ensembl_gene_id"],
        "previous_symbols": record["prev_symbol"],
        "alias_symbols": record["alias_symbol"],
    })
    return annotation


def main():
    if not INPUT_FILE.is_file():
        raise FileNotFoundError(
            f"Missing file: {INPUT_FILE}\nRun PrepareRAGData.py first."
        )

    results = pd.read_csv(
        INPUT_FILE, dtype={"gene": str}, keep_default_na=False
    )

    if results.empty:
        raise ValueError("The retained gene table is empty.")

    if not {"gene", "comparison_id"}.issubset(results.columns):
        raise ValueError("Missing gene or comparison_id column.")

    if results["gene"].str.strip().eq("").any():
        raise ValueError("Empty gene labels found.")

    reference = load_reference()
    current, historical, records = build_indexes(reference)

    # Annotate each distinct label once, even across several cancers.
    annotations = pd.DataFrame([
        annotate_gene(gene, current, historical, records)
        for gene in sorted(results["gene"].unique())
    ])

    enriched = results.merge(
        annotations,
        on="gene",
        how="left",
        validate="many_to_one",
    )

    # Different original labels may map to the same official gene.
    # Keep these rows separate and flag them instead of merging statistics.
    enriched["mapping_collision"] = False
    mapped = enriched["mapping_status"].eq("mapped")

    enriched.loc[mapped, "mapping_collision"] = (
        enriched.loc[mapped]
        .duplicated(["comparison_id", "hgnc_id"], keep=False)
    )

    enriched["needs_review"] = (
        ~mapped | enriched["mapping_collision"]
    )

    annotations.to_csv(DATA_DIR / "gene_annotations.csv", index=False)
    enriched.to_csv(
        DATA_DIR / "retained_genes_annotated.csv", index=False
    )
    enriched.loc[enriched["needs_review"]].to_csv(
        DATA_DIR / "genes_to_review.csv", index=False
    )

    print("\nMapping summary — distinct original gene labels:")
    print(annotations["mapping_status"].value_counts().to_string())

    print(
        "\nCancer–gene rows requiring review:",
        int(enriched["needs_review"].sum()),
    )
    print(f"\nOutputs saved in: {DATA_DIR}")


if __name__ == "__main__":
    main()