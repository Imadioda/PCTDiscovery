# Differential Expression Analysis with limma

This module compares primary tumor tissue with anatomically compatible normal tissue from the same TCGA diagnostic cohort. Python scripts annotate and prepare the samples; R performs the statistical analysis. Precomputed results are provided for inspection without rerunning the workflow.

**Current status:** all 11 exported comparisons completed successfully, with 16,335 expression features tested in each comparison. Results remain exploratory. Completion of the calculations is not equivalent to biological or clinical validation.

## Folder layout

This README belongs in the repository's `resultats_DE_preparation/` folder.

| Relative path | Contents |
|---|---|
| `preparation/GDCClassification.py` | GDC specimen metadata retrieval and anatomical annotation |
| `preparation/PrepareDEComparisons.py` | Eligibility filtering and per-comparison exports |
| `resultats_DE/` | The complete Python export folder, renamed for this repository layout |
| `AnalyseLimma.R` | Local R analysis entry point |
| `Limma_Colab.ipynb` | Colab workflow with an embedded R script |
| `runlimma/` | Files and comparison subfolders from the completed reference run |

`resultats_DE/` must contain `plan_analyses.csv` directly, alongside the comparison folders such as `LUAD_lung/`. Preserve all relative paths in the plan. Each comparison contains `expression.tsv.gz`, `samples.csv`, `comparison.json`, and `genes_excluded.csv`.

For `runlimma/`, copy the **contents** of the completed `run_...` directory here, so that `runlimma/synthese_analyses.csv` exists directly. Preserve `sessionInfo.txt`, `execution.txt`, `plan_entree.csv`, all comparison folders and any warnings. The parent session's `console.log` and executed R script may also be included for provenance. Do not mix files from different runs.

## Three ways to use this module

### 1. Inspect the published results

Start with [the run summary](runlimma/synthese_analyses.csv), then open the corresponding cancer folder under `runlimma/`. No Python or R installation is needed to inspect CSV tables and PDF plots.

| Output per comparison | Meaning |
|---|---|
| `resultats_complets.csv` | All tested genes, model effects, TREAT p-values and BH-adjusted p-values |
| `genes_retenus.csv` | Genes selected at the configured FDR threshold |
| `genes_surexprimes.csv` | Selected genes with positive tumor-minus-normal effects |
| `genes_sous_exprimes.csv` | Selected genes with negative tumor-minus-normal effects |
| `echantillons_utilises.csv` | Retained samples and patient/group metadata |
| `matrice_design.csv` | Statistical design matrix |
| `qc_echantillons.csv` | Per-sample expression summaries |
| `qc_PCA_coordonnees.csv`, `qc_PCA.pdf` | PCA coordinates and visualization |
| `qc_distributions.pdf` | Sample expression distributions |
| `qc_moyenne_variance.pdf` | Mean–variance diagnostic |
| `volcano.pdf`, `MA_plot.pdf` | Effect size and statistical-result visualizations |
| `parametres_et_limites.txt` | Analysis settings and limitations |

### 2. Rerun limma from the prepared data

Install R, then run the following in an **R console**:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install("limma", ask = FALSE, update = FALSE)
install.packages(
    c("statmod", "data.table", "R.utils", "matrixStats"),
    repos = "https://cloud.r-project.org"
)
```

The original R script needs `limma` and `statmod`; the optimized Colab variant additionally uses `data.table`, `R.utils`, and `matrixStats`. Installing this full set supports both variants. Consult the published `sessionInfo.txt` for the reference run's actual package versions; these installation commands do not lock versions.

From the **repository root**, run in a terminal:

```bash
Rscript resultats_DE_preparation/AnalyseLimma.R \
  resultats_DE_preparation/resultats_DE \
  resultats_DE_preparation/reruns
```

The first argument is the directory containing `plan_analyses.csv`. The second is the output base directory. The script creates a timestamped run inside that output directory. Explicit arguments override the original local path defaults and keep reruns separate from the published `runlimma/` reference.

**Do not unzip `expression.tsv.gz`.** R reads these compressed matrices directly. Neither `dataset_enrichi.pkl` nor `dataset_DE_filtre.pkl` is used by R. No classifier checkpoint is required.

### 3. Reproduce the Python preparation

Install the preparation dependencies:

```bash
python -m pip install pandas numpy
```

Start from the project's cleaned `dtf_modele` DataFrame, saved as `dtf_modele.pkl` at the repository root. It must contain the original Xena log-expression values, before model-specific standardization, with sample identifiers, `sample_type`, `_primary_disease`, and numeric expression columns. A sample-named index is also accepted.

Run from the repository root:

```bash
python resultats_DE_preparation/preparation/GDCClassification.py \
  --input dtf_modele.pkl \
  --output resultats_DE_preparation/resultats_GDC

python resultats_DE_preparation/preparation/PrepareDEComparisons.py \
  --input resultats_DE_preparation/resultats_GDC/dataset_enrichi.pkl \
  --output resultats_DE_preparation/resultats_DE_rebuilt
```

Use a new output directory for the second command: preparation rejects nonempty output folders. To analyze these rebuilt data, use `resultats_DE_rebuilt` as the R input path instead of `resultats_DE`.

The first script saves the **enriched** DataFrame, `dataset_enrichi.pkl`, and its accompanying `genes.json`. Keep these together. The second saves the **filtered** union, `dataset_DE_filtre.pkl`, and the per-comparison exports. The pickle files are Python intermediates, not portable R inputs. Only load trusted pickle files.

The second script does **not** automatically create a ZIP archive.

## Optional Google Colab execution

Open `Limma_Colab.ipynb` in Google Colab, using a standard Python/CPU runtime. It installs R and runs the embedded R script. A GPU is not used by this workflow.

To create the input ZIP, run from the repository root:

```bash
cd resultats_DE_preparation
zip -r donnees_limma.zip resultats_DE -x '*.pkl' '*/.DS_Store'
```

Run notebook cells in order, enable Google Drive saving if desired, and upload this archive when prompted. The notebook discovers the directory containing `plan_analyses.csv`, so the export folder can be called `resultats_DE`. It reads inputs from the runtime's local disk and can save outputs to Drive.

The notebook uses its **embedded script**, not the adjacent `AnalyseLimma.R`. Changes to the local R file do not automatically update the notebook.

**There is no automatic resume in the current notebook.** Rerunning the analysis cell starts all comparisons again in a new session directory. The download cell packages only the directory currently held in `SESSION_DIR`, not all historical sessions. Previously saved Drive directories are preserved. An interrupted comparison must be recomputed unless a separate checkpoint/restart workflow is implemented.

The reference execution took approximately ten hours on the user's runtime; this is not a runtime guarantee. `duplicateCorrelation` was the main observed bottleneck. Progress messages report elapsed time, not an estimated completion time. Cloud resources and session lifetime are not guaranteed.

## Cohort construction

The annotation script queries the GDC API and specimen-linked BCR Biotab supplements. It records sample-level anatomical metadata and source provenance. It does not substitute a patient's cancer diagnosis or case-level primary site for a missing specimen location. Ambiguous or conflicting annotations are retained for review or excluded from the analysis exports.

Comparisons are predefined, not all possible cancer–organ combinations. Each requires:

- The expected TCGA project and disease label.
- `Primary Tumor` or `Solid Tissue Normal` sample type.
- Documented compatible anatomy for the corresponding group.
- At least 10 independent patients per group after filtering.
- No repeated eligible observation for the same patient within the same group.

A normal and tumor sample from the same patient can both remain. Patient IDs are exported to account for this dependence in R. Normal controls come from the same diagnostic cohort, not necessarily from healthy donors. Different cancer subtypes in one organ remain separate comparisons.

Of 17 predefined comparisons, 11 were exported and completed in the reference run. Inspect `resultats_DE/plan_analyses.csv` and `decisions_echantillons.csv` for exclusions and skipped comparisons.

## Statistical method

The input is already normalized and on the **`log2(norm_value + 1)`** scale. This workflow does not use DESeq2, voom, another log transformation, imputation, or StandardScaler.

For each comparison, the script:

1. Checks expression/metadata alignment, identifiers, finite values, group counts and gene counts.
2. Produces expression distributions and PCA using up to 500 most variable genes for visualization only.
3. Fits a linear model with normal and tumor groups and a **Tumor − Normal** contrast.
4. Uses `duplicateCorrelation` and a patient-block model when paired samples exist. Fewer than three pairs trigger design review rather than silently assuming independence. With no pairs, an independent-patient model is used.
5. Applies robust, trend-aware empirical Bayes **TREAT**, testing whether the absolute effect exceeds **1 unit on the input log-expression scale**.
6. Adjusts p-values with Benjamini–Hochberg and selects genes at **FDR ≤ 0.05**, separately within each comparison.

The effect is a model-estimated difference on the input scale, not necessarily the log2 ratio of raw counts. With patient correlation accounted for, it can differ from the simple difference between the exported unweighted group means. All exported genes are tested; the PCA subset is not a gene-selection step for limma.

## Reference results

Counts below come from the completed `synthese_analyses.csv`. All statuses are `completed_exploratory`; all comparisons tested 16,335 genes.

| Comparison | Normal | Tumor | Paired patients | Upregulated | Downregulated |
|---|---:|---:|---:|---:|---:|
| THCA_thyroid | 58 | 494 | 58 | 501 | 639 |
| BRCA_breast | 101 | 1,062 | 96 | 736 | 1,653 |
| PRAD_prostate | 43 | 484 | 43 | 115 | 356 |
| KIRC_kidney | 32 | 511 | 32 | 1,267 | 1,223 |
| KIRP_kidney | 30 | 289 | 30 | 510 | 1,362 |
| LUAD_lung | 58 | 495 | 51 | 857 | 1,354 |
| LUSC_lung | 46 | 455 | 44 | 1,317 | 2,308 |
| LIHC_liver | 43 | 370 | 43 | 519 | 705 |
| COAD_colon | 41 | 432 | 37 | 681 | 1,718 |
| STAD_stomach | 28 | 415 | 25 | 501 | 511 |
| UCEC_endometrium | 17 | 482 | 5 | 772 | 1,342 |

Genes can be selected in multiple comparisons; summing these counts does not give the number of unique genes across cancers.

## Quality review and limitations

A detailed output-consistency review was performed for **LUAD**, including identifier checks, counts, recomputation of BH-adjusted p-values, selection/direction consistency, and visual QC. The input expression matrix was not independently refitted during that review. This is not a complete validation of all 11 comparisons.

The tumor sample `TCGA-86-A4P8-01` lies near normal samples on the first two LUAD PCs. It is flagged for metadata review, not automatically removed or relabelled. The first two PCs capture 24.4% of variance in the 500-gene visualization subset.

Remaining limits include unmodelled age, sex, additional batch and tumor-purity covariates; possible tissue-composition/subsite differences; inherited upstream Xena processing; and exclusion of genes absent from the common 16,335-feature matrix. The minimum group size is a workflow threshold, not a power calculation. FDR correction is within comparisons, not across the entire cancer collection.

Differential expression establishes association, not causality, predictive usefulness, or a validated therapeutic target. Biological annotation, pathway enrichment, sensitivity analyses and literature-supported interpretation remain future work.

## Sources

- [UCSC Xena expression dataset](https://xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net)
- [GDC API documentation](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/)
- [TCGA sample type codes](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)
- [limma](https://bioconductor.org/packages/limma/)
- [TREAT methodology](https://pmc.ncbi.nlm.nih.gov/articles/PMC2654802/)
- [Google Colab FAQ](https://research.google.com/colaboratory/faq.html)

Keep source URLs, GDC annotation manifests, eligibility decisions and `sessionInfo.txt` with the corresponding data/run. If matrices are distributed outside Git, document their download location before publishing this README.
