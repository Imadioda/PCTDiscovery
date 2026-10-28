# Phase 1 schema — version 0.1.0

`graph.json` separates normalized nodes, graph assertions and per-source evidence. IDs are strings, never row positions. All exported object IDs are stable for fixed inputs and source snapshots. Retrieval/build timestamps are excluded from the core graph except retained source provenance, so rebuilding the same snapshots is deterministic.

## Nodes

| Type | Namespace | Interpretation |
|---|---|---|
| Gene | `NCBIGene:` or unresolved-crosswalk `HGNC:` | Normalized gene locus; a PubTator protein mention does not create a separate protein entity |
| Disease | `MESH:` | Normalized disease concept; broad neoplasms remain broad |
| Chemical | `MESH:` | Chemical concept, not necessarily a drug or an approved therapy |
| Species | `NCBITaxon:` | Species mention |
| CellLine | `Cellosaurus:CVCL_` | Cell-line identity |
| Variant | `dbSNP:rs` | dbSNP-linked variant; other forms are currently quarantined |
| Publication | `PMID:` | Source publication |
| StudyComparison | `PTCComparison:` | Our cohort/contrast, deliberately distinct from a disease ontology concept |

Each node has `id`, `type`, `namespace`, `labels`, `xrefs`, `attributes`. Labels are sorted alternatives, not evidence that similarly named nodes are identical. Available taxon features can be recorded through HGNC cross-references; this phase does not produce learned embeddings.

## Edges

Each edge has `id`, `subject`, `predicate`, `object`, `status`, `directed` and `evidence_ids`. IDs are hashes of the assertion identity, not scores. Multiple source assertions support one graph edge where identities agree; their evidence records remain separate.

- `directed=false`: a relation treated as symmetric; endpoint IDs are sorted deterministically.
- `directed=null`: biological direction unknown; original source roles remain in evidence. Do not interpret storage order as causality.
- `directed=true`: an explicitly directed dataset observation, such as gene→study comparison.

Predicates `pubtator:*` preserve extraction semantics such as association or positive/negative correlation. They are not mapped to stronger labels such as "causes" or "validated target". `ptc:differentially_expressed_in` carries expression results only.

Statuses: `automatically_extracted`, `dataset_observation`, `documented`, `predicted`. Only the first two are currently populated. Review status lives at evidence level; extraction and human validation are not conflated.

## Evidence

Extraction evidence retains PMID, original PubTator relation ID and object, source role order, extraction scope, score if present, snapshot URL/hash/time and supporting annotation IDs/passages where available. The original API can identify a supporting mention pair through `nodes[].role`; its `refid` is not always an annotation ID. Both the current PubTator3 layout and conventional BioC annotation-reference layouts are handled.

Limma evidence retains original gene label, HGNC ID, table row, effect, FDR, cohort details and scale `log2(normalized_value + 1)`. The build manifest records the CSV hash. Do not convert these effects into exact raw-count fold changes.

Corpus comparison/gene labels are retrieval provenance only. Tissue, species, negation, experimental intervention, assay type, isoforms and causal direction are not automatically inferred from article membership. Where unresolved, preserve the original text and require review.

## Review and future extensions

`issues.csv` quarantines normalization failures; it is not a list of biologically unimportant entities. `review_queue.csv` supports source-level human evaluation. Graph integrity validation and scientific validation are separate gates.

A future curated-resource importer can add UniProt proteins, gene–protein cross-references, reactions, pathways and curated interactions with resource version and source evidence. A later prediction exporter must write status `predicted` with model/checkpoint/split identity and scores; it must never masquerade as literature evidence.
