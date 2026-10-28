# PTCDiscovery — Knowledge Graph, phase 1

A separate, evidence-aware preparation module for the planned GNN and GraphRAG stages. It reads the existing RAG data snapshot, retrieves PubTator 3 annotations/relations, normalizes stable identifiers and exports a heterogeneous graph with a record for each source assertion.

**Status: working phase-1 pilot, not a scientifically validated knowledge base.** No GNN is trained and the existing RAG application is not modified.

## Quick start

Place this `knowledge_graph` directory at the root of PTCDiscovery. From it, use **Python 3.11 or newer**:

```bash
cd knowledge_graph
python -m unittest discover -s tests -v
```

The module uses only the Python standard library: no API key, GPU, database server, PyTorch or extra pip dependencies are required at this stage. Use your usual interpreter (`python3` instead of `python` if needed).

### Pilot with existing RAG data

Locate the actual `rag_data` directory extracted from your version 1 Release. It must contain `retained_genes_annotated.csv` and `literature_all_cancers/articles.json` and `gene_article_links.csv`. The original HGNC snapshot at `references/hgnc_complete_set.txt` is strongly recommended for gene-ID reconciliation.

Replace `/absolute/path/to/rag_data` below with your own path. Do not copy old scripts into one folder for this module: it already uses package-relative imports and explicit input paths.

```bash
python -m ptc_kg run \
  --rag-data /absolute/path/to/rag_data \
  --comparison LUAD_lung \
  --pmids 36405247 37350390 38309204 \
  --limit 3 \
  --collection outputs/luad_collection \
  --output outputs/luad_graph
```

This is the real-data pilot tested during development. It retrieves those three publications and imports **all eligible retained LUAD expression rows**, not only genes in the three papers. Explicit PMIDs override automatic corpus selection, while preserving any known corpus context. Numerical ordering of PMIDs defines the capped selection; it is not a relevance or quality ranking.

Without `--pmids`, the input corpus's eligible publications are selected, filtered by comparison and capped by `--limit` (default 10). Use `--comparison all` for all comparisons once the pilot is reviewed. It is not necessary to start with a large collection.

### Enrich the collection through PubTator search

```bash
python -m ptc_kg run \
  --query 'osimertinib lung adenocarcinoma' \
  --pages 1 --limit 2 \
  --collection outputs/chemical_collection \
  --output outputs/chemical_graph
```

This publication-only example needs no RAG files. A query can also be combined with existing seeds; their PMID union is sorted numerically before applying the limit. Search is bounded to 1–5 pages. A comparison label is a selection filter/context, not an automatic disease assignment to every result. PubTator indexing/search results may evolve.

### Separate collection and construction

Use `collect` with the collection options to create a `collection_manifest.json`, then:

```bash
python -m ptc_kg build \
  --collection outputs/luad_collection \
  --rag-data /absolute/path/to/rag_data \
  --output outputs/luad_graph_rebuild

python -m ptc_kg validate outputs/luad_graph_rebuild/graph.json
```

Use a **new or empty graph output directory** for each build. Successful cached responses are reused. To replay the same `run` command without network access, add `--offline` and choose a new graph output directory. Every requested response must already be cached. `build` itself makes no network request.

Failed requests are recorded; `collect` exits nonzero. Construction refuses incomplete collections unless you explicitly pass `--allow-partial`. Retry collection to resume from the successful cache rather than deleting it. Raw input hashes must match when the same RAG snapshot is supplied at build time.

## Files and responsibilities

| File | Role |
|---|---|
| `ptc_kg/pubtator.py` | Rate-limited REST access, bounded retries, content-hashed JSON snapshots and offline replay |
| `ptc_kg/normalize.py` | Namespace-aware normalization and conservative HGNC–NCBI Gene crosswalk |
| `ptc_kg/graph.py` | Node identity, separate assertions/evidence, expression observations, integrity validation |
| `ptc_kg/pipeline.py` | Existing-corpus selection, exclusions, collection/build manifests and exports |
| `ptc_kg/common.py` | Strict JSON, CSV, booleans, hashes and atomic file replacement |
| `ptc_kg/__main__.py` | `collect`, `build`, `run`, `validate` commands |
| `ptc_kg/exclusions.json` | Version 1 association-specific exclusion policy, packaged with the code |
| `tests/test_phase1.py` | Offline tests for identity, extraction, provenance, failures and deterministic rebuilds |
| `SOURCE_AUDIT.md` | Existing-repository inspection and reuse decisions |
| `SCHEMA.md` | Node/edge/evidence meanings and interpretation rules |

## Outputs

The collection directory contains `cache/` raw response snapshots, `collection_progress.json` and `collection_manifest.json`. Collection progress is saved after each publication. Interrupted runs retain earlier successful responses.

The graph directory contains:

- `graph.json`: authoritative full graph including raw relation objects, supporting passages, entity mentions and issues.
- `nodes.csv`: normalized entities, node types, aliases, cross-references and available attributes.
- `edges.csv`: typed, evidence-backed assertions, grouped by endpoints/predicate/status.
- `evidence.csv`: source-level records, extraction scores and original role order; nested values are JSON cells.
- `mentions.csv`: annotated text spans, locations and verification status.
- `issues.csv`: unresolved/ambiguous identifiers, mapping flags and other review items.
- `review_queue.csv`: one row per extracted edge–evidence record, with blank reviewer/label/notes columns.
- `validation_report.json`: integrity checks, graph sizes, type counts and issue counts.
- `build_manifest.json`: source/input/code hashes, version, environment and limitations.

The normalizer does not guess identities from gene symbols. Valid PubTator3 Cellosaurus IDs such as `CVCL:1561` map to `Cellosaurus:CVCL_1561`. Unsupported variant identifiers remain review items instead of being silently discarded without a trace.

## Evidence rules

1. PubTator/BioREx relations are **automatically extracted and unreviewed**, even when their source is a published article. Their extraction scores are not calibrated biological probabilities.
2. `dataset_observation` edges store our limma results against **study-comparison nodes**. They do not establish gene–disease causation, molecular interaction or target validity.
3. `documented` and `predicted` are reserved separate statuses. This phase does not populate them or promote PubTator extractions to curated evidence.
4. A PMID's presence in a cancer-specific corpus does not mean every relation in that paper is cancer-specific. Corpus provenance is stored separately from biological context, which remains `not_resolved`.
5. Co-mentions are not converted to interactions. Missing edges are not negative labels.
6. PubTator's gene/protein mentions generally normalize to **NCBI Gene**. This phase does not infer UniProt proteins, isoforms, gene→protein edges, pathways or approved treatments. Dedicated resources will be required for those entity/relation types.
7. Relation role order is preserved, but no causal direction is inferred from it. Known symmetric labels are represented symmetrically; other relation directions remain unknown until reviewed.
8. Whole source passages can be retained as context. They are not claimed to be sentence-level proof of every relation. Original BioC locations/objects remain available for inspection.
9. HGNC crosswalks must be unique in both directions. Missing or ambiguous mappings are flagged. No human species identity is inferred merely because a paper mentions a human cancer; human taxonomy attributes are added only through the supplied approved HGNC crosswalk.
10. Retraction/expression-of-concern exclusions reuse the supplied article metadata. New publications and old metadata still need current scientific/retraction review before modeling.

## Validation gate before the GNN

The automated gate checks IDs, referential integrity, evidence attachment, numeric values and provenance. It does **not** establish scientific correctness. Review `review_queue.csv` against each source and classify examples as supported, unsupported, ambiguous or context-mismatched. Keep the article, raw relation and context with each review. This first version exports the queue; it does not yet ingest reviews to alter graph status.

Before training, define one prediction task and its eligible evidence tier. Deduplicate the same biological fact across papers, group inverse/duplicate relations when splitting, and consider temporal or independent-source holdouts. Review text features and limma-derived attributes for leakage into the target labels. Do not train on an automatically extracted relation while also using its full assertion or equivalent inverse as an input feature. Unobserved relations need an explicit negative-sampling policy.

A graph schema with multiple types does not imply enough data for a GNN. The small pilot is intended to validate ingestion. Structural QC, biological review and a task-specific split come before R-GCN/GAT/GraphSAGE comparisons; GSL and GraphRAG remain later phases.

## Sources

- PubTator3 API: https://www.ncbi.nlm.nih.gov/research/pubtator3/api
- BioC export: https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson
- PubTator paper: https://doi.org/10.1093/nar/gkae235
- NCBI API example repository: https://github.com/ncbi-nlp/pubtator-gpt
- HGNC reference: the original version 1 `hgnc_complete_set.txt` snapshot.

Raw downloaded passages remain third-party publication content. Preserve provenance and review applicable reuse terms before redistributing a full corpus. No private API credentials are needed by this module.

## Evidence review

See [REVIEW.md](REVIEW.md) for the evidence-level HTML/CSV review workflow and `review-export` / `review-import` commands. Reviews preserve the original graph and do not automatically create training labels.
