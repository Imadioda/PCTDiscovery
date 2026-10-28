# Validation record — 21 October 2026

## Automated tests

Command, executed from `knowledge_graph/`:

```bash
python -m unittest discover -s tests -v
```

**25 tests passed**, offline and without third-party Python dependencies. Coverage includes:

- supported BioC wrappers and both relation endpoint formats;
- identifier namespaces, ambiguous IDs, normalized-ID conflicts and crosswalk collisions;
- gene/protein distinction and native Cellosaurus normalization;
- no interaction from co-mentions, no causal direction from source role order;
- distinct PMID evidence, within-document relation deduplication and review flags;
- original mention span verification and invalid extraction scores;
- broken graph references/evidence detection;
- source-corpus exclusions and retracted-publication filtering;
- bounded network retries, wrong-PMID responses, cache hash checks and offline misses;
- expression observation separation and deterministic graph exports.

## Real NCBI API tests

Both BioC export and PubTator search were called successfully. The tests below are ingestion/integrity checks, not biological validation or measures of PubTator precision.

| Check | LUAD + existing results | Query-enrichment pilot |
|---|---:|---:|
| Requested/collected publications | 3 / 3 | 2 / 2 |
| Failed publication requests | 0 | 0 |
| Nodes | 2,192 | 15 |
| Automatically extracted graph edges | 21 | 13 |
| Limma dataset-observation edges | 2,166 | 0 |
| Total source evidence records | 2,187 | 17 |
| Verified/retained mention records | 117 | 78 |
| Review issues | 46 | 12 |
| Integrity report | valid | valid |

LUAD PMIDs: **36405247, 37350390, 38309204**. RAG inputs came from repository Release `1.0`. Of its retained LUAD rows, 45 mapping-review rows were quarantined; one accepted HGNC gene lacked a unique usable NCBI Gene crosswalk and remained an HGNC node with a review issue.

Enrichment query: **`osimertinib lung adenocarcinoma`**, one search page, cap 2 after numeric PMID ordering. Exported PMIDs: **33410885, 33919291**. The second pilot includes Chemical nodes. Its issues are five ambiguous annotation IDs and seven unresolved relation endpoints, all retained for review rather than guessed. Relation types include association, cotreatment and correlations; none is promoted to validated treatment efficacy.

The two graph counts are separate test outputs, not a merged validated database. Graph edges may aggregate distinct source records, so evidence and edge counts can differ. All extracted relations remain unreviewed.

## Offline reproducibility

The real LUAD snapshot was rebuilt with the same input files and `--offline`. `graph.json` was byte-identical to its online-snapshot rebuild. Raw retrieval timestamps belong to the stored source snapshot, while each build has its own manifest timestamp.

Source snapshots can be replayed locally once collected. Raw API caches and complete article passages are not bundled in this source archive. The supplied small validation tables are inspection results; they are not a self-contained corpus for replay. Re-run the documented pilot with your Release data to fetch those snapshots.

## Evidence available for review

`validation/luad/` and `validation/chemical/` contain the validation reports and the source-level review queue from each pilot. They include IDs and scores, not full publication texts. The source article and original raw export must be consulted before filling a human review label.

## Remaining gate

No GNN, GSL or GraphRAG has been implemented or evaluated here. No comparison of model architectures is claimed. Before proceeding, review a representative sample of graph relations, establish context-specific label rules and define leakage-resistant splits. This phase provides the ingestion foundation and an auditable review queue.
