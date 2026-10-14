# Differential Expression Analysis with limma

This module compares primary tumor tissue with anatomically compatible normal tissue from the same TCGA diagnostic cohort. Python scripts annotate and prepare the samples; R performs the statistical analysis. Precomputed results are provided for inspection without rerunning the workflow.

**Current status:** all 11 exported comparisons completed successfully, with 16,335 expression features tested in each comparison. Results remain exploratory. Completion of the calculations is not equivalent to biological or clinical validation.

## Repository files and Release archives

This README belongs in `differential_expression/`.

| Repository path (relative to this folder) | Contents |
|---|---|
| `preparation/GDCClassification.py` | GDC metadata retrieval and anatomical annotation |
| `preparation/PrepareDEComparisons.py` | Filtering and per-comparison exports |
| `AnalyseLimma.R` | Local R analysis |
| `Limma_Colab.ipynb` | Colab workflow with an embedded R script |
| [plan_analyses.csv](plan_analyses.csv) | Preparation eligibility decisions and comparison counts |
| [synthese_analyses.csv](synthese_analyses.csv) | Completed-run summary |

The two large archives are distributed as assets in this repository's **Releases** section, not as ordinary Git files. Open the release containing the differential-expression data and results and download the attached assets, not just GitHub's automatically generated source-code archive.

- **Prepared-data archive:** the original Python export directory, including `plan_analyses.csv` and each comparison's `expression.tsv.gz`, `samples.csv`, `comparison.json`, and `genes_excluded.csv`.
- **Results archive:** the original completed run/session directory, including per-cancer result tables and plots, `synthese_analyses.csv`, run parameters and software information.

**Archive contents retain their original names and nesting. Do not flatten or rename their internal folders.** Extract each archive to its own local directory. The input directory for R is whichever extracted directory directly contains `plan_analyses.csv`. The reference result directory is whichever directly contains `synthese_analyses.csv`; it may be nested inside `session_.../run_...`.

The two small CSV copies in this repository are for quick consultation. The standalone `plan_analyses.csv` is not enough to run R: its referenced expression and metadata files must also be present in their original relative locations.

## Three ways to use this module

### 1. Inspect the published results

Start with [the run summary](synthese_analyses.csv), then download and extract the results archive from Releases and open the corresponding cancer folder inside the completed run. No Python or R installation is needed to inspect CSV tables and PDF plots.

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
Rscript differential_expression/AnalyseLimma.R \
  "/absolute/path/to/extracted_preparation_folder" \
  differential_expression/reruns
```

Replace the example first argument with the actual **directory** containing `plan_analyses.csv`, not the CSV file itself. The second is the output base directory. The script creates a timestamped run inside that output directory. Explicit arguments override the original local path defaults and keep reruns separate from the extracted reference results.

**Do not unzip `expression.tsv.gz`.** R reads these compressed matrices directly. Neither `dataset_enrichi.pkl` nor `dataset_DE_filtre.pkl` is used by R. No classifier checkpoint is required.

### 3. Reproduce the Python preparation

Install the preparation dependencies:

```bash
python -m pip install pandas numpy
```

Start from the project's cleaned `dtf_modele` DataFrame, saved as `dtf_modele.pkl` at the repository root. It must contain the original Xena log-expression values, before model-specific standardization, with sample identifiers, `sample_type`, `_primary_disease`, and numeric expression columns. A sample-named index is also accepted.

Run from the repository root:

```bash
python differential_expression/preparation/GDCClassification.py \
  --input dtf_modele.pkl \
  --output differential_expression/resultats_GDC

python differential_expression/preparation/PrepareDEComparisons.py \
  --input differential_expression/resultats_GDC/dataset_enrichi.pkl \
  --output differential_expression/resultats_DE_rebuilt
```

Use a new output directory for the second command: preparation rejects nonempty output folders. To analyze these rebuilt data, use `resultats_DE_rebuilt` as the R input path instead of the downloaded input directory.

The first script saves the **enriched** DataFrame, `dataset_enrichi.pkl`, and its accompanying `genes.json`. Keep these together. The second saves the **filtered** union, `dataset_DE_filtre.pkl`, and the per-comparison exports. The pickle files are Python intermediates, not portable R inputs. Only load trusted pickle files.

The second script does **not** automatically create a ZIP archive.

## Optional Google Colab execution

Open `Limma_Colab.ipynb` in Google Colab, using a standard Python/CPU runtime. It installs R and runs the embedded R script. A GPU is not used by this workflow.

Use the **prepared-data ZIP from the Release** as the notebook input. Run cells in order, enable Drive saving if desired, and upload that archive when prompted. Do not upload the results archive. The notebook recursively locates `plan_analyses.csv`, so original wrapping directories can remain intact. The uploaded ZIP must contain exactly one preparation plan.

It reads input matrices from the runtime's local disk and can save outputs to Drive. If the Release's input asset is already a ZIP, no manual recompression is needed.

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

Of 17 predefined comparisons, 11 were exported and completed in the reference run. Inspect [plan_analyses.csv](plan_analyses.csv) and the extracted `decisions_echantillons.csv` for exclusions and skipped comparisons.

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

Keep source URLs, GDC annotation manifests, eligibility decisions and `sessionInfo.txt` with the corresponding data/run. The large data and result archives are available in the repository’s Releases section; the scripts and summary CSVs remain in this folder.
