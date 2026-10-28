# Evidence review — phase 1 validation

This step reviews whether a source supports an extracted relation. It does not establish causality, therapeutic efficacy, or training readiness. Existing graph files are immutable; reviews are separate, versioned-by-input snapshots.

## Run

From `knowledge_graph/` (Python 3.11+, no external dependencies):

```bash
python -m unittest discover -s tests -v
python -m ptc_kg review-export --graph outputs/luad_graph/graph.json --output outputs/luad_review
```

Open `outputs/luad_review/review.html` in a browser. Edit the corresponding `decisions.csv` in a CSV editor, preserving UTF-8, comma delimiters, column names and identifiers. The HTML is read-only and contains the available context, source metadata, extraction roles and PubMed links. Empty context means the source needs to be consulted directly; do not guess.

For each evidence row:

- `supported`: the source supports the exact normalized entities and relation. Supply reviewer, notes, a short supporting quotation, a source locator (e.g. abstract sentence or full-text section), and biological context (species, tissue/cancer, experimental setting where available).
- `unsupported`: the extraction is not supported, is misnormalized, or misrepresents the source. Supply reviewer and explanatory notes. This is NOT a negative biological link for training.
- `uncertain`: the source or context is insufficient or ambiguous. Supply reviewer and notes.
- Leave all review fields blank when not reviewed. Never delete rows.

Source locator and quotation are reviewer assertions; the importer checks completeness, not their truth. Check negation, direction, species, cell-line versus patient results, and association versus mechanism. A positive extraction review alone never changes `automatically_extracted` to `documented`.

```bash
python -m ptc_kg review-import --graph outputs/luad_graph/graph.json --decisions outputs/luad_review/decisions.csv --output outputs/luad_reviewed
```

Outputs: `reviews.json`, `edge_review_summary.csv`, `review_report.json`. Imports reject a mismatched graph hash, edited identity columns, unknown or duplicated evidence, missing rows, and incomplete decisions. Outputs must be new/empty directories. Source graph is never overwritten. Multiple reviewers should use separate review bundles; adjudication is not implemented yet. Mixed outcomes across publications are flagged for inspection, not automatically interpreted as biological contradiction.

## Included real pilot bundles

`review_pilots/luad/`: 21 evidence records, 21 extracted edges.
`review_pilots/chemical/`: 17 evidence records, 13 extracted edges.

Each contains `graph.json`, `review.html`, blank `decisions.csv`, `review_manifest.json`. These graphs are the original collected pilot snapshots, not new evidence. Their raw API caches are not bundled. LUAD also includes the existing 2,166 expression observations, which are not part of extraction review. Import an edited pilot CSV against the graph in that same folder.

No scientific decisions have been prefilled. This is a small non-random development sample, not an estimate of extraction accuracy across cancers. Broader validation needs a prespecified sampling plan, relation/cancer coverage, independent review and adjudication. A future GNN also needs a task definition, a defensible positive-label policy, negative sampling that does not equate unknown with false, and leakage-resistant splits. `training_ready` therefore remains false.

## Verification

35 automated tests pass (25 phase-1 tests plus 10 review tests). Real pilot export produced 21 and 17 review records respectively. Tests use synthetic decisions; no test decision is a scientific annotation.
