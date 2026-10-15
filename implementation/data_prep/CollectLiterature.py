from pathlib import Path
from datetime import datetime, timezone
from hashlib import sha256
from html import unescape
import json
import re
import time

from html.parser import HTMLParser


import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "rag_data" / "retained_genes_annotated.csv"
OUTPUT_DIR = BASE_DIR / "rag_data" / "literature_LUAD"
CACHE_DIR = OUTPUT_DIR / "cache"

COMPARISON = "LUAD_lung"
GENES_PER_DIRECTION = 10
ARTICLES_PER_GENE = 5
CANDIDATES_PER_QUERY = 50

API_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# Keep the first collection specific to lung adenocarcinoma.
CANCER_TERMS = [
    "lung adenocarcinoma",
    "pulmonary adenocarcinoma",
    "adenocarcinoma of the lung",
    "LUAD",
]


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def save_json(path, data):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(
        json.dumps(data, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    temporary.replace(path)


def save_csv(path, table):
    temporary = path.with_suffix(".tmp")
    table.to_csv(temporary, index=False)
    temporary.replace(path)


def read_boolean(series):
    values = series.astype(str).str.strip().str.lower()
    if not values.isin(["true", "false"]).all():
        raise ValueError(f"Invalid boolean values in {series.name}")
    return values.eq("true")


class AbstractTextParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.parts = []

    def handle_data(self, data):
        self.parts.append(data)

    def handle_starttag(self, tag, attrs):
        self.parts.append(" ")

    def handle_endtag(self, tag):
        self.parts.append(" ")


def clean_text(text):
    parser = AbstractTextParser()
    parser.feed(text or "")
    parser.close()
    return " ".join("".join(parser.parts).split())


def text_query(term):
    term = term.replace('"', " ").strip()
    return f'(TITLE:"{term}" OR ABSTRACT:"{term}")'


def create_query(row):
    # No automatic expansion to potentially ambiguous aliases.
    gene_terms = [row["approved_symbol"]]
    if row["gene_name"].strip():
        gene_terms.append(row["gene_name"])

    gene_query = " OR ".join(text_query(term) for term in gene_terms)
    cancer_query = " OR ".join(text_query(term) for term in CANCER_TERMS)

    return (
        f"({gene_query}) AND ({cancer_query}) "
        "AND HAS_ABSTRACT:Y AND SRC:MED"
    )


def create_session():
    session = requests.Session()

    retries = Retry(
        total=4,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        respect_retry_after_header=True,
    )

    session.mount("https://", HTTPAdapter(max_retries=retries))
    session.headers.update({
        "User-Agent": "CancerExpressionLiteraturePilot/0.1"
    })
    return session


def search_articles(session, query):
    parameters = {
        "query": query,
        "format": "json",
        "resultType": "core",
        "pageSize": CANDIDATES_PER_QUERY,
    }

    key = sha256(
        json.dumps(parameters, sort_keys=True).encode()
    ).hexdigest()
    cache_file = CACHE_DIR / f"{key}.json"

    if cache_file.exists():
        return json.loads(cache_file.read_text(encoding="utf-8"))

    response = session.get(
        API_URL,
        params=parameters,
        timeout=(15, 90),
    )
    response.raise_for_status()
    data = response.json()

    if "hitCount" not in data or "resultList" not in data:
        raise ValueError("Unexpected Europe PMC response.")

    payload = {
        "retrieved_at_utc": utc_now(),
        "request_url": response.url,
        "parameters": parameters,
        "response": data,
    }
    save_json(cache_file, payload)

    # Modest request frequency for this pilot.
    time.sleep(0.5)
    return payload


def select_genes():
    if not INPUT_FILE.is_file():
        raise FileNotFoundError(f"Missing input: {INPUT_FILE}")

    table = pd.read_csv(INPUT_FILE, keep_default_na=False)

    required = {
        "comparison_id", "gene", "approved_symbol", "gene_name",
        "hgnc_id", "mapping_method", "mapping_status", "needs_review",
        "selected", "effect_log2_scale", "fdr_bh_treat",
    }
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    eligible = table.loc[
        table["comparison_id"].eq(COMPARISON)
        & table["mapping_status"].eq("mapped")
        & table["mapping_method"].eq("approved_symbol")
        & ~read_boolean(table["needs_review"])
        & read_boolean(table["selected"])
    ].copy()

    for column in ["effect_log2_scale", "fdr_bh_treat"]:
        eligible[column] = pd.to_numeric(
            eligible[column], errors="raise"
        )

    if eligible["hgnc_id"].eq("").any():
        raise ValueError("Missing HGNC identifiers.")

    if eligible["hgnc_id"].duplicated().any():
        raise ValueError("Duplicate HGNC identifiers in this comparison.")

    eligible["absolute_effect"] = eligible["effect_log2_scale"].abs()
    eligible = eligible.sort_values(
        ["absolute_effect", "fdr_bh_treat", "approved_symbol"],
        ascending=[False, True, True],
    )

    up = eligible.loc[
        eligible["effect_log2_scale"] > 0
    ].head(GENES_PER_DIRECTION)

    down = eligible.loc[
        eligible["effect_log2_scale"] < 0
    ].head(GENES_PER_DIRECTION)

    if len(up) < GENES_PER_DIRECTION or len(down) < GENES_PER_DIRECTION:
        raise ValueError(
            "Not enough eligible genes for 10 up / 10 down."
        )

    return pd.concat([up, down], ignore_index=True)

def main():
    global COMPARISON, CANCER_TERMS, OUTPUT_DIR, CACHE_DIR

    # Specific search terms for each comparison.
    # Avoid broad terms that would mix different cancer subtypes.
    comparisons = {
        "LUAD_lung": [
            "lung adenocarcinoma",
            "pulmonary adenocarcinoma",
            "adenocarcinoma of the lung",
            "LUAD",
        ],
        "LUSC_lung": [
            "lung squamous cell carcinoma",
            "pulmonary squamous cell carcinoma",
            "squamous cell carcinoma of the lung",
        ],
        "BRCA_breast": [
            "invasive breast carcinoma",
            "invasive breast cancer",
            "breast carcinoma",
            "breast cancer",
        ],
        "THCA_thyroid": [
            "thyroid carcinoma",
            "thyroid cancer",
            "papillary thyroid carcinoma",
        ],
        "PRAD_prostate": [
            "prostate adenocarcinoma",
            "prostatic adenocarcinoma",
            "prostate cancer",
        ],
        "KIRC_kidney": [
            "clear cell renal cell carcinoma",
            "renal clear cell carcinoma",
            "ccRCC",
        ],
        "KIRP_kidney": [
            "papillary renal cell carcinoma",
            "papillary renal carcinoma",
        ],
        "LIHC_liver": [
            "hepatocellular carcinoma",
            "hepatocellular cancer",
        ],
        "COAD_colon": [
            "colon adenocarcinoma",
            "colonic adenocarcinoma",
            "colon cancer",
        ],
        "STAD_stomach": [
            "gastric adenocarcinoma",
            "stomach adenocarcinoma",
            "gastric cancer",
        ],
        "UCEC_endometrium": [
            "endometrial carcinoma",
            "endometrial cancer",
            "endometrioid endometrial carcinoma",
            "uterine corpus endometrial carcinoma",
        ],
    }

    if not INPUT_FILE.is_file():
        raise FileNotFoundError(f"Missing input: {INPUT_FILE}")

    table = pd.read_csv(INPUT_FILE, keep_default_na=False)

    required = {
        "comparison_id", "gene", "approved_symbol", "gene_name",
        "hgnc_id", "mapping_method", "mapping_status", "needs_review",
        "selected", "effect_log2_scale", "fdr_bh_treat",
    }
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    available = set(table["comparison_id"])
    missing_comparisons = set(comparisons) - available
    if missing_comparisons:
        raise ValueError(
            "Missing comparisons in the input: "
            + ", ".join(sorted(missing_comparisons))
        )

    eligible = table.loc[
        table["mapping_status"].eq("mapped")
        & table["mapping_method"].eq("approved_symbol")
        & ~read_boolean(table["needs_review"])
        & read_boolean(table["selected"])
    ].copy()

    for column in ["effect_log2_scale", "fdr_bh_treat"]:
        eligible[column] = pd.to_numeric(
            eligible[column], errors="raise"
        )

    effects = eligible["effect_log2_scale"]
    if not effects.map(
        lambda value: float("-inf") < value < float("inf")
    ).all():
        raise ValueError("Non-finite expression effects.")

    if not eligible["fdr_bh_treat"].between(0, 1).all():
        raise ValueError("Invalid FDR values.")

    if eligible["hgnc_id"].eq("").any():
        raise ValueError("Missing HGNC identifiers.")

    if eligible.duplicated(["comparison_id", "hgnc_id"]).any():
        raise ValueError("Duplicate cancer–HGNC pairs.")

    # Build the selection before making API requests.
    selections = []
    selection_summary = []

    for comparison in comparisons:
        group = eligible.loc[
            eligible["comparison_id"].eq(comparison)
        ].copy()

        group["absolute_effect"] = group["effect_log2_scale"].abs()
        group = group.sort_values(
            ["absolute_effect", "fdr_bh_treat", "approved_symbol"],
            ascending=[False, True, True],
        )

        up = group.loc[
            group["effect_log2_scale"] > 0
        ].head(GENES_PER_DIRECTION)

        down = group.loc[
            group["effect_log2_scale"] < 0
        ].head(GENES_PER_DIRECTION)

        selections.extend([up, down])
        selection_summary.append({
            "comparison_id": comparison,
            "selected_up": len(up),
            "selected_down": len(down),
        })

    selected = pd.concat(selections, ignore_index=True)
    if selected.empty:
        raise ValueError("No eligible genes.")

    # Separate outputs from the original LUAD pilot.
    OUTPUT_DIR = BASE_DIR / "rag_data" / "literature_all_cancers"
    CACHE_DIR = OUTPUT_DIR / "cache"
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # Reuse existing pilot responses without modifying the original cache.
    pilot_cache = BASE_DIR / "rag_data" / "literature_LUAD" / "cache"

    if pilot_cache.is_dir():
        import shutil

        copied = 0
        for source in pilot_cache.glob("*.json"):
            target = CACHE_DIR / source.name
            if not target.exists():
                shutil.copy2(source, target)
                copied += 1
        print(f"Pilot cache responses imported: {copied}")

    save_csv(OUTPUT_DIR / "selected_genes.csv", selected)
    save_csv(
        OUTPUT_DIR / "selection_summary.csv",
        pd.DataFrame(selection_summary),
    )

    save_json(OUTPUT_DIR / "collection_manifest.json", {
        "started_at_utc": utc_now(),
        "input_sha256": sha256(INPUT_FILE.read_bytes()).hexdigest(),
        "comparison_queries": comparisons,
        "genes_per_direction": GENES_PER_DIRECTION,
        "articles_per_gene": ARTICLES_PER_GENE,
        "candidates_per_query": CANDIDATES_PER_QUERY,
        "selection": "Largest absolute effects; exact approved symbols",
        "article_order": "Europe PMC default relevance",
        "scope": "PubMed-indexed records with abstracts",
        "review_status": "Candidates only; relevance not validated",
        "limitations": [
            "Search terms do not guarantee matching cancer histology.",
            "Gene mentions do not establish functional evidence.",
            "Only the first candidate page is inspected per query.",
            "Missing results do not establish absence of literature.",
        ],
    })

    documents = {}
    links = []
    summary = []

    link_columns = [
        "comparison_id", "gene", "approved_symbol", "hgnc_id",
        "article_id", "search_rank", "query", "retrieved_at_utc",
        "review_status",
    ]

    def save_progress():
        save_json(
            OUTPUT_DIR / "articles.json",
            list(documents.values()),
        )
        save_csv(
            OUTPUT_DIR / "gene_article_links.csv",
            pd.DataFrame(links, columns=link_columns),
        )
        save_csv(
            OUTPUT_DIR / "collection_summary.csv",
            pd.DataFrame(summary),
        )

    print("\nSelection:")
    print(pd.DataFrame(selection_summary).to_string(index=False))
    print(f"\nTotal cancer–gene pairs: {len(selected)}")

    with create_session() as session:
        for position, row in enumerate(
            selected.to_dict("records"), start=1
        ):
            COMPARISON = row["comparison_id"]
            CANCER_TERMS = comparisons[COMPARISON]
            symbol = row["approved_symbol"]

            print(
                f"[{position}/{len(selected)}] "
                f"{COMPARISON} — {symbol}",
                flush=True,
            )

            query = create_query(row)

            try:
                payload = search_articles(session, query)
                data = payload["response"]
                candidates = data["resultList"].get("result", [])

                count = 0
                seen = set()

                for rank, article in enumerate(candidates, start=1):
                    abstract = clean_text(article.get("abstractText"))
                    source = article.get("source", "")
                    identifier = article.get("id", "")

                    if not abstract or not identifier or not source:
                        continue

                    article_id = f"{source}:{identifier}"
                    if article_id in seen:
                        continue
                    seen.add(article_id)

                    publication_types = (
                        (article.get("pubTypeList") or {})
                        .get("pubType", [])
                    )
                    if isinstance(publication_types, str):
                        publication_types = [publication_types]

                    flags = [
                        value for value in publication_types
                        if "retract" in value.lower()
                        or "expression of concern" in value.lower()
                    ]

                    # Store one document per Europe PMC identifier.
                    if article_id not in documents:
                        documents[article_id] = {
                            "article_id": article_id,
                            "source": source,
                            "source_id": identifier,
                            "pmid": article.get("pmid", ""),
                            "pmcid": article.get("pmcid", ""),
                            "doi": article.get("doi", ""),
                            "title": clean_text(article.get("title")),
                            "abstract": abstract,
                            "authors": article.get("authorString", ""),
                            "year": article.get("pubYear", ""),
                            "publication_types": publication_types,
                            "publication_flags": flags,
                            "url": (
                                f"https://europepmc.org/article/"
                                f"{source}/{identifier}"
                            ),
                            "retrieved_at_utc": payload["retrieved_at_utc"],
                            "content_scope": "abstract",
                            "review_status": "pending",
                        }

                    # Preserve each cancer–gene–article association.
                    links.append({
                        "comparison_id": COMPARISON,
                        "gene": row["gene"],
                        "approved_symbol": symbol,
                        "hgnc_id": row["hgnc_id"],
                        "article_id": article_id,
                        "search_rank": rank,
                        "query": query,
                        "retrieved_at_utc": payload["retrieved_at_utc"],
                        "review_status": "pending",
                    })

                    count += 1
                    if count >= ARTICLES_PER_GENE:
                        break

                summary.append({
                    "comparison_id": COMPARISON,
                    "gene": row["gene"],
                    "approved_symbol": symbol,
                    "hits": int(data["hitCount"]),
                    "articles_collected": count,
                    "status": "completed" if count else "no_articles",
                    "error": "",
                })

                print(
                    f"  {count} collected / "
                    f"{data['hitCount']} search hits",
                    flush=True,
                )

            except (
                requests.RequestException,
                ValueError,
                KeyError,
            ) as error:
                summary.append({
                    "comparison_id": COMPARISON,
                    "gene": row["gene"],
                    "approved_symbol": symbol,
                    "hits": "",
                    "articles_collected": 0,
                    "status": "error",
                    "error": str(error),
                })
                print(f"  Error: {error}", flush=True)

            save_progress()

    report = pd.DataFrame(summary)
    report["has_articles"] = report["articles_collected"].gt(0)
    report["failed"] = report["status"].eq("error")

    coverage = (
        report.groupby("comparison_id", as_index=False)
        .agg(
            genes_searched=("gene", "size"),
            genes_with_articles=("has_articles", "sum"),
            gene_article_links=("articles_collected", "sum"),
            failed_queries=("failed", "sum"),
        )
    )

    save_csv(OUTPUT_DIR / "coverage_by_cancer.csv", coverage)

    print("\nCoverage:")
    print(coverage.to_string(index=False))
    print(f"\nUnique articles: {len(documents)}")
    print(f"Gene–article links: {len(links)}")
    print(f"Outputs: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()