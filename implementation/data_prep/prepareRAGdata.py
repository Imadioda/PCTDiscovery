from pathlib import Path
import pandas as pd


# Folder containing the eleven cancer result folders.
LIMMA_DIR = Path(
    "/run_limma"
)

OUTPUT_DIR = Path(__file__).resolve().parent / "rag_data"


def prepare_rag_data():
    if not LIMMA_DIR.is_dir():
        raise FileNotFoundError(f"Folder not found: {LIMMA_DIR}")

    files = sorted(LIMMA_DIR.rglob("resultats_complets.csv"))

    if not files:
        raise FileNotFoundError(
            f"No resultats_complets.csv found in {LIMMA_DIR}"
        )

    required = {
        "gene",
        "comparison_id",
        "cancer",
        "effect_log2_scale",
        "fdr_bh_treat",
        "selected",
        "direction",
    }

    tables = []

    for file in files:
        table = pd.read_csv(file)

        missing = required - set(table.columns)
        if missing:
            raise ValueError(f"{file}: missing columns {sorted(missing)}")

        if table.empty or table[list(required)].isna().any().any():
            raise ValueError(f"{file}: empty table or missing required values")

        # Read booleans explicitly: bool("False") would incorrectly be True.
        selected = table["selected"].astype(str).str.lower().str.strip()
        if not selected.isin(["true", "false"]).all():
            raise ValueError(f"{file}: invalid selected values")

        table["selected"] = selected.eq("true")

        for column in ["effect_log2_scale", "fdr_bh_treat"]:
            table[column] = pd.to_numeric(table[column], errors="raise")

        effect = table["effect_log2_scale"]
        if not effect.map(
            lambda value: float("-inf") < value < float("inf")
        ).all():
            raise ValueError(f"{file}: non-finite expression effects")

        if not table["fdr_bh_treat"].between(0, 1).all():
            raise ValueError(f"{file}: invalid FDR values")

        if table["comparison_id"].nunique() != 1:
            raise ValueError(f"{file}: expected one comparison per file")

        # Preserve the original gene labels and all statistical columns.
        table["source_file"] = file.relative_to(LIMMA_DIR).as_posix()
        tables.append(table)

    results = pd.concat(tables, ignore_index=True)

    # Prevent accidentally combining several runs of the same analysis.
    duplicated = results.duplicated(["comparison_id", "gene"], keep=False)
    if duplicated.any():
        raise ValueError(
            "Duplicate cancer–gene entries found. "
            "Choose a single limma run folder."
        )

    retained = results.loc[results["selected"]].copy()
    retained["absolute_effect"] = retained["effect_log2_scale"].abs()

    retained = retained.sort_values(
        ["comparison_id", "absolute_effect", "fdr_bh_treat", "gene"],
        ascending=[True, False, True, True],
    )

    retained["effect_rank"] = (
        retained.groupby("comparison_id").cumcount() + 1
    )

    summary = (
        results.groupby(["comparison_id", "cancer"], as_index=False)
        .agg(
            genes_tested=("gene", "size"),
            genes_retained=("selected", "sum"),
        )
    )

    # A new folder prevents accidental replacement of an existing export.
    OUTPUT_DIR.mkdir(parents=True, exist_ok=False)

    results.to_csv(OUTPUT_DIR / "all_results.csv", index=False)
    retained.to_csv(OUTPUT_DIR / "retained_genes.csv", index=False)
    summary.to_csv(OUTPUT_DIR / "cancer_summary.csv", index=False)

    print(summary.to_string(index=False))
    print(f"\nComparisons: {len(summary)}")
    print(f"Retained cancer–gene pairs: {len(retained)}")
    print(f"Unique retained gene labels: {retained['gene'].nunique()}")
    print(f"Saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    prepare_rag_data()