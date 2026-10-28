# Existing components inspected

Repository: https://github.com/Imadioda/PCTDiscovery

Baseline commit examined: `ea01a5c1762712ad8352c8d7322168ac40f865bc`.

The repository's actual folders are lowercase `implementation/data_prep/` and `implementation/answers_build/`; the new module uses neither the older flat-import assumption nor capitalization from earlier planning.

| Existing component | Reuse decision |
|---|---|
| `implementation/data_prep/AnnotateGenes.py` | Reuse the released annotation CSV and HGNC reference snapshot, including `needs_review`/mapping collision flags. Do not rerun or bypass this mapping review. |
| `implementation/data_prep/CollectLiterature.py` | Reuse `articles.json`, publication flags and the cancer–gene–article table. Its links are collection candidates, not validated biological interactions. |
| `implementation/answers_build/LiteratureSearch.py` | Preserve the specific LUAD–CLDN18 article exclusion and retraction/concern filtering as explicit new-module policy. No importing the v1 script at runtime. |
| `implementation/answers_build/BuildRAGContext.py` | Reuse statistical field names; do not use generated contexts or LLM answers as biological labels. |
| `AskRAG.py`, `LLMProviders.py`, `AppRAG.py` | No modification. Phase 1 does not call a paid LLM. |
| Existing model and differential-expression scripts | No modification or execution. Use released statistical results only. |

Data snapshot inspected: Release `1.0`, asset `rag_data.zip`. Only the article corpus, association CSV, annotated retained results and HGNC snapshot were extracted for testing. Their per-file hashes are recorded in the generated manifests. The asset URL is `https://github.com/Imadioda/PCTDiscovery/releases/download/1.0/rag_data.zip`.

The API was tested against real NCBI responses before writing the parser. Current exports use a `PubTator3` root and document-level relations with normalized `infons.role1`/`role2` objects. The parser also handles conventional BioC collections/references without inventing endpoints for unsupported structures.

All existing tracked files remain unchanged. The deliverable is an added `knowledge_graph/` module; no upstream push or deployment is performed.
