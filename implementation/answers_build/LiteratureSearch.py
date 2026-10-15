from pathlib import Path
import json

import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity


BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "rag_data" / "literature_all_cancers"

# Exclusion specific to this cancer–gene association:
# the abstract describes gastric CLDN18 findings, not LUAD evidence.
EXCLUDED_LINKS = {
    ("LUAD_lung", "CLDN18", "MED:42562391"),
}


class LiteratureSearch:
    def __init__(self, data_dir=DATA_DIR):
        data_dir = Path(data_dir)

        with (data_dir / "articles.json").open(encoding="utf-8") as handle:
            articles = json.load(handle)

        links = pd.read_csv(
            data_dir / "gene_article_links.csv",
            dtype=str,
            keep_default_na=False,
        )

        required = {
            "article_id", "comparison_id", "approved_symbol",
            "hgnc_id", "review_status",
        }
        if not required.issubset(links.columns):
            raise ValueError("Missing required link columns.")

        identifiers = [article["article_id"] for article in articles]
        if len(identifiers) != len(set(identifiers)):
            raise ValueError("Duplicate article identifiers.")

        if links.duplicated(
            ["comparison_id", "hgnc_id", "article_id"]
        ).any():
            raise ValueError("Duplicate cancer–gene–article links.")

        missing = set(links["article_id"]) - set(identifiers)
        if missing:
            raise ValueError(f"Missing articles: {sorted(missing)}")

        self.cancers = sorted(links["comparison_id"].unique())

        excluded_ids = set()

        for article in articles:
            labels = (
                article.get("publication_flags", [])
                + article.get("publication_types", [])
            )
            if any(
                "retract" in label.lower()
                or "expression of concern" in label.lower()
                for label in labels
            ):
                excluded_ids.add(article["article_id"])

        keep_link = [
            (
                row.comparison_id,
                row.approved_symbol,
                row.article_id,
            ) not in EXCLUDED_LINKS
            for row in links.itertuples(index=False)
        ]

        self.links = links.loc[
            keep_link
        ].loc[
            lambda frame: ~frame["article_id"].isin(excluded_ids)
        ].copy()

        linked_ids = set(self.links["article_id"])
        self.documents = [
            article for article in articles
            if article["article_id"] in linked_ids
            and article.get("abstract", "").strip()
        ]

        if not self.documents:
            raise ValueError("No documents available for indexing.")

        # One document per abstract: preserve its surrounding context.
        texts = [
            f"{article.get('title', '')}\n{article['abstract']}"
            for article in self.documents
        ]

        self.vectorizer = TfidfVectorizer(
            lowercase=True,
            stop_words="english",
            ngram_range=(1, 2),
            sublinear_tf=True,
            max_features=60000,
        )
        self.matrix = self.vectorizer.fit_transform(texts)

        print(f"Indexed articles: {len(self.documents)}")
        print(f"Excluded flagged documents: {len(excluded_ids)}")
        print(f"Cancers: {len(self.cancers)}")

    def search(self, question, cancer, gene=None, top_k=5):
        if cancer not in self.cancers:
            raise ValueError(
                f"Unknown cancer: {cancer}. Options: {self.cancers}"
            )
        if not question.strip():
            raise ValueError("The question cannot be empty.")
        if top_k < 1:
            raise ValueError("top_k must be positive.")

        # Filter associations BEFORE ranking documents.
        associations = self.links.loc[
            self.links["comparison_id"].eq(cancer)
        ]

        if gene:
            associations = associations.loc[
                associations["approved_symbol"]
                .str.casefold()
                .eq(gene.strip().casefold())
            ]

        allowed_ids = set(associations["article_id"])

        positions = [
            index for index, article in enumerate(self.documents)
            if article["article_id"] in allowed_ids
        ]

        if not positions:
            return []

        query_vector = self.vectorizer.transform([question])

        if query_vector.nnz == 0:
            return []

        scores = cosine_similarity(
            query_vector,
            self.matrix[positions],
        ).ravel()

        ranking = sorted(
            range(len(positions)),
            key=lambda index: (
                -float(scores[index]),
                self.documents[positions[index]]["article_id"],
            ),
        )

        results = []

        for local_index in ranking[:top_k]:
            score = float(scores[local_index])
            if score <= 0:
                continue

            article = self.documents[positions[local_index]]
            article_links = associations.loc[
                associations["article_id"].eq(article["article_id"])
            ]

            results.append({
                "article_id": article["article_id"],
                "comparison_id": cancer,
                "linked_genes": sorted(
                    article_links["approved_symbol"].unique().tolist()
                ),
                "retrieval_score": round(score, 4),
                "title": article.get("title", ""),
                "abstract": article["abstract"],
                "year": article.get("year", ""),
                "doi": article.get("doi", ""),
                "url": article.get("url", ""),
                "association_review_status": sorted(
                    article_links["review_status"].unique().tolist()
                ),
                "limitation": (
                    "Candidate literature association; "
                    "retrieval score is not scientific evidence."
                ),
            })

        return results


if __name__ == "__main__":
    engine = LiteratureSearch()

    results = engine.search(
        question="GPM6A proliferation migration PI3K AKT",
        cancer="LUAD_lung",
        gene="GPM6A",
        top_k=3,
    )

    for position, result in enumerate(results, start=1):
        print(f"\n{position}. {result['title']}")
        print(f"Search score: {result['retrieval_score']}")
        print(f"Genes: {', '.join(result['linked_genes'])}")
        print(f"Source: {result['url']}")
        print(result["abstract"])

    if not results:
        print("No matching documents found.")

    output_file = DATA_DIR / "search_example.json"
    output_file.write_text(
        json.dumps(results, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"\nExample saved to: {output_file}")