# PCTDiscovery

**Version 0 · Active development · Independent bioinformatics research project**

A research workflow combining TCGA gene-expression classification, anatomically matched tumor–normal differential expression, and planned literature-supported biomarker and therapeutic-target exploration.

## Overview

I am developing this project as a first-year bioinformatics master's student to gain experience across public-data acquisition, quality control, machine learning, statistical analysis and biological interpretation.

The current implementation includes:

1. A logistic regression model distinguishing normal tissue from primary tumor tissue.
2. A PyTorch multilayer perceptron ranking 32 cancer types for samples predicted as tumors.
3. GDC specimen-level anatomical annotation and preparation of tumor–normal comparisons.
4. Completed limma differential-expression analyses for 11 cancer–tissue comparisons, with downloadable tables and QC plots.

The inputs are **bulk RNA-seq tissue samples**, not individual cells. A normal-tissue prediction does not imply that the donor is cancer-free. This is an exploratory research prototype, not a clinically validated diagnostic tool.

The project will be updated regularly. Literature retrieval, biological interpretation and candidate-target prioritization are still under development.

## Inspiration

This project is inspired by the research objectives of the Sanofi–Owkin collaboration: using AI in oncology for biomarker identification, therapeutic-target discovery and treatment-response research.

It is an independent implementation using public data. **It is not affiliated with or endorsed by Sanofi or Owkin**, and does not reproduce their proprietary systems, federated-learning infrastructure or research results.

Reference: [Sanofi–Owkin announcement, November 18, 2021](https://www.sanofi.com/en/media-room/press-releases/2021/2021-11-18-06-30-00-2336966).

## Data and preprocessing

Expression and phenotype data originate from the **TCGA Pan-Cancer Atlas**, distributed through UCSC Xena.

| Input | Filename |
|---|---|
| Expression | `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz` |
| Phenotypes | `TCGA_phenotype_denseDataOnlyDownload.tsv` |

Sources: [expression dataset](https://xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net) and [Pan-Cancer Atlas hub](https://xenabrowser.net/datapages/?host=https%3A%2F%2Fpancanatlas.xenahubs.net).

Expression values are already normalized and on the **`log2(norm_value + 1)`** scale. No additional log transformation is applied. After transposition and merging, rows represent samples, with expression features in columns and `sample`, `sample_type`, `_primary_disease` as identifiers/annotations. Disease labels describe the diagnostic cohort and may also appear on its normal tissue samples.

| Processing stage | Samples | Expression features |
|---|---:|---:|
| Merged matrix after duplicate-name handling | 11,060 | 20,531 |
| After exploratory zero-variance filtering | 11,060 | 20,312 |
| Complete-feature matrix | 11,060 | 16,335 |
| Normal/primary-tumor modeling subset | 10,438 | 16,335 |
| Primary tumors for cancer classification | 9,701 | 16,335 |

The complete prepared table has 16,338 columns including the three identifier/annotation columns. Expression features are referred to as genes for convenience. Duplicate-name processing retained `SLC35E2_1` and `SLC35E2_2`; their mapping/provenance still needs explicit verification.

### Missing-data decision

After the initial variance filtering, **3,588,487 missing values** remained, distributed across only five exact missingness profiles:

| Profile | Samples | Missing genes per sample |
|---|---:|---:|
| 0 | 9,310 | 0 |
| 1 | 33 | 1,224 |
| 2 | 908 | 1,236 |
| 3 | 211 | 3,547 |
| 4 | 598 | 2,805 |

Missingness was strongly associated with cancer cohorts. Removing every incomplete sample would eliminate the acute myeloid leukemia cohort and almost all ovarian, stomach and esophageal samples.

Version 0 therefore removes **3,977 incomplete features** and retains **16,335 complete features without imputation**. This preserves samples at the cleaning stage; model-specific sample-type restrictions are applied afterward. Excluded genes are not assumed biologically unimportant. The technical origin of the missingness blocks remains unresolved.

## Prediction pipeline

### Binary model

The binary dataset contains **9,701 primary tumors and 737 solid-tissue normal samples**. Other sample types are outside the current prediction scope. This follows the annotations: lymphomas labelled `Primary Tumor` can be included, while blood-derived leukemia samples are excluded.

- scikit-learn pipeline with training-only variance filtering and standardization.
- L2-regularized logistic regression, L-BFGS, balanced class weights.
- Validation selection among `C = 0.01, 0.1, 1.0`; selected `C = 1.0`.
- Decision threshold 0.5, fixed before the internal test evaluation.
- Saved preprocessing/model artifact in `.joblib` format.

### Cancer classifier

A PyTorch **MLP**, not an RNN, is trained on true primary tumors in the training partition, independently of binary predictions.

Architecture: **16,335 → 128 → 64 → 32**.

- Hidden layers use BatchNorm, ReLU and dropout (0.3 and 0.2).
- A separate StandardScaler is fitted on network training samples only.
- Unweighted cross-entropy on logits; softmax at inference.
- AdamW, learning rate 0.001, weight decay 0.0001, batch size 64.
- Early stopping based on validation loss, patience 10 epochs.
- Best checkpoint at epoch 9; training stopped at epoch 19.
- Weights, preprocessing parameters, feature order and class mapping saved in `.pth` format.

At inference, the second model runs only on samples predicted as tumors. It returns the top three cancer labels and softmax scores. Scores are not yet calibrated and are not established clinical probabilities. The cancer scores are not multiplied by the binary score.

## Evaluation

Both models reuse patient-grouped partitions. All samples from one patient remain in one partition. The split attempts to preserve normal/tumor proportions; cancer-class coverage is checked afterward.

| Partition | Normal | Tumor | Total |
|---|---:|---:|---:|
| Train | 442 | 5,820 | 6,262 |
| Validation | 147 | 1,940 | 2,087 |
| Internal test | 148 | 1,941 | 2,089 |

Validation results: binary balanced accuracy **98.30%**, ROC-AUC **0.9876**; cancer-classifier accuracy **95.88%**, balanced accuracy **92.75%**, macro-F1 **92.73%**, top-3 accuracy **99.64%**.

### Final internal test

These results use all 2,089 reserved test samples, not only the randomly displayed examples.

| Evaluation | Metric | Result |
|---|---|---:|
| Binary | Accuracy | 99.33% |
| Binary | Balanced accuracy | 98.39% |
| Binary | ROC-AUC | 0.9994 |
| Binary | Missed tumors | 10 / 1,941 |
| Binary | False positives among normals | 4 / 148 |
| Complete chain, true tumors | Detected and correctly identified, top-1 | 1,863 / 1,941 (95.98%) |
| Complete chain, true tumors | Detected with true cancer in top-3 | 1,926 / 1,941 (99.23%) |
| Complete chain, all samples | Accuracy | 96.07% |
| Complete chain, all samples | Macro-F1 across 33 labels | 93.69% |

Binary confusion matrix:

| Actual / Predicted | Normal | Tumor |
|---|---:|---:|
| Normal | 144 | 4 |
| Tumor | 10 | 1,931 |

A tumor rejected at the first stage counts as an error in both complete-chain tumor metrics. Overall final accuracy is 2,007 / 2,089. Rectal adenocarcinoma remains a major weakness, with 8 / 28 correctly classified. Small class sizes also limit interpretation: cholangiocarcinoma has only three test samples.

### Evaluation limitations

- Complete-feature selection and initial variance filtering used the full dataset before splitting. Training-only model preprocessing does not eliminate this earlier information leakage; the full evaluation remains exploratory.
- Upstream normalization and cohort/tissue/batch effects may influence performance.
- This is internal testing, not external validation. The inspected test set must not become a tuning set while still being described as untouched.
- Calibration, uncertainty handling and external generalization remain unassessed.

## Differential expression: completed analysis module

See **[the differential-expression README](differential_expression/README.md)** for commands, dependencies, input/output formats, per-cancer counts and limitations.

The workflow annotates specimen anatomy through GDC, then prepares predefined cancer–tissue comparisons. It retains compatible primary tumors and solid-tissue normal controls from the same diagnostic cohort, requiring at least ten patients per group. Paired normal/tumor samples are preserved and modeled through patient blocks.

**Eleven comparisons completed:** THCA, BRCA, PRAD, KIRC, KIRP, LUAD, LUSC, LIHC, COAD, STAD and UCEC. These are 11 cancer comparisons, not 11 distinct organs. Each tested 16,335 expression features.

The analysis uses **limma with robust, trend-aware TREAT**, a minimum effect of 1 on the existing log-expression scale, and **BH FDR ≤ 0.05 within each comparison**. No DESeq2, voom, additional log transformation or model-standardized values are used. PCA is a QC visualization, not the differential-expression test.

The run summary records successful completion for all comparisons. LUAD also underwent detailed output-consistency and visual QC checks, confirming 857 selected upregulated and 1,354 selected downregulated features. This review did not independently refit the expression matrix or validate all other cancers. A tumor sample near normal samples in the LUAD PCA is documented for follow-up.

Results remain exploratory: age, sex, additional batch and tumor purity are not included as covariates in this version. Anatomical matching cannot guarantee identical cell composition. Differential expression does not establish causality or a therapeutic target.

This analysis branch is separate from predictive evaluation. If its gene lists are later used as classifier features, selection must be repeated within training data rather than using all-sample differential-expression results to tune the previously tested models.

## Repository organization and execution

The differential-expression code and small summary tables are under `differential_expression/`; the large prepared-data and completed-run archives are attached to the repository's **Release**.

| Repository path | Purpose |
|---|---|
| `differential_expression/README.md` | Detailed instructions |
| `differential_expression/preparation/` | Two Python preparation scripts |
| `differential_expression/AnalyseLimma.R` | Local R execution |
| `differential_expression/Limma_Colab.ipynb` | Optional cloud execution |
| [differential_expression/plan_analyses.csv](differential_expression/plan_analyses.csv) | Eligibility and export summary |
| [differential_expression/synthese_analyses.csv](differential_expression/synthese_analyses.csv) | Completed-run summary |

Open this repository's **Releases** section to download the two attached archives. One contains the prepared expression matrices and metadata; the other contains the completed limma results. Their original internal folders and filenames are unchanged. GitHub's automatically generated source-code ZIP is not a substitute for these separately attached data assets.

After extracting the prepared-data archive, pass R the directory directly containing `plan_analyses.csv`. After extracting the results archive, find `synthese_analyses.csv` inside the reference run, possibly under `session_.../run_...`. Do not flatten the archive structure. See the module README for commands.

For the original classification scripts, the documented local workflow uses a flat working directory containing:

- `SampleDataClean.py`, `DataMerge.py`, `ExpressionDataClean.py` for preprocessing.
- `LogModelTrain.py` for binary training.
- `NeuralNtwr.py` for the network definition and `NrnModelTrain.py` for training.
- `PipelineModel.py` for inference and `PipelineTest.py` for final evaluation.

Some earlier copies use `ModeleBinaire.py` or `LancerPrediction.py`; use the corresponding entry point in your checkout. If these files are grouped into presentation folders such as `processing/`, `models/` and `final_test/`, the current flat-layout copies must be placed together for execution unless their imports and paths have been updated for packages. Merely moving files does not update imports or data/model paths. Training and evaluation artifact locations must match the paths configured in the scripts.

### Python classification setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The development environment used Python 3.13 with Miniconda on macOS. Core classification dependencies include NumPy, pandas, scikit-learn, PyTorch and joblib. The dependency list is not yet a fully locked environment.

Download the Xena input files and place them at the paths expected by the preprocessing scripts. For the flat local layout, use the project root. Prepare `dtf_modele`, train the binary model to save partitions and preprocessing, then train the neural network using those same partitions. Run `PipelineTest.py` as a normal Python script, not through pytest, once the models and threshold are fixed.

Loading a classifier requires the matching architecture, feature order, scaler, class mapping and compatible dependencies. Unless model checkpoints are distributed separately, train them before running the evaluation. `NeuralNtwr.py` should define the network without launching preprocessing or training at import time.

### Differential-expression setup

To inspect results, consult `differential_expression/synthese_analyses.csv` or download the results archive from Releases; no calculation is needed. To rerun limma, use the prepared exports and the module's R instructions. To reproduce cohort preparation, run its Python scripts first.

The scripts do not require a model checkpoint for differential expression. Python pickle intermediates are not R inputs; the R script reads `plan_analyses.csv`, `expression.tsv.gz` and `samples.csv`.

## Runtime and distribution

Large expression matrices consume substantial RAM. A dense 11,060 × 20,531 float64 matrix alone is approximately 1.82 GB before DataFrame overhead and copies. Preprocessing imports may repeat expensive loading in the current classification implementation.

The reference differential-expression run took approximately ten hours on the user's runtime. Patient-correlation estimation was the main observed bottleneck. Actual duration depends on hardware, sample counts and package versions. The current limma code uses CPU resources; selecting a GPU does not accelerate it.

The supplied Colab notebook saves each run in a new session directory and optionally writes to Google Drive. It has **no automatic resume**. Restarting the analysis cell begins all comparisons again; its download cell packages only the currently selected session. See the module README before rerunning a completed or interrupted analysis.

Keep raw Xena files, caches, large pickle files, temporary reruns and model checkpoints out of routine Git history. The prepared-data and completed-results archives are distributed through Releases with their original structure preserved. The preparation plan and run summary are directly available under `differential_expression/`. Ensure `.gitignore` does not exclude these two small published tables.

## Progress and next steps

- [x] Merge expression and phenotype data and investigate structured missingness.
- [x] Build a complete 16,335-feature matrix without imputation.
- [x] Train and save a binary logistic regression and 32-class PyTorch MLP.
- [x] Link both models and evaluate the complete chain on the internal test set.
- [x] Retrieve GDC specimen metadata and annotate anatomical provenance.
- [x] Prepare predefined tumor–normal comparisons with explicit eligibility rules.
- [x] Complete 11 limma analyses and retain results, QC plots and run metadata.
- [x] Review output consistency and visual QC in detail for LUAD.
- [ ] Review the remaining comparisons, flagged samples and potential confounding.
- [ ] Pin dependencies and verify feature identifiers, including duplicate-name handling.
- [ ] Improve execution paths, caching and checkpoint/resume support.
- [ ] Rebuild predictive evaluation with training-only feature selection and independent assessment.
- [ ] Compare simpler multiclass baselines and investigate class-specific errors.
- [ ] Assess calibration and external generalization.
- [ ] Annotate selected genes and perform pathway/functional enrichment analysis.
- [ ] Explore molecular subgroups and stability of predictive and differential signatures.
- [ ] Prioritize candidate genes using statistical, biological and network evidence.
- [ ] Build a RAG layer retrieving literature on genes, proteins, pathways and documented interventions.
- [ ] Generate source-backed reports and develop an exploration interface.

The long-term objective is to connect computed results with documented biological evidence. Candidate therapeutic targets will remain research hypotheses requiring further validation. The current repository does not yet implement the RAG or target-prioritization stages.

## Further references

- [GDC API](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/)
- [limma package and documentation](https://bioconductor.org/packages/limma/)
- [TREAT methodology](https://pmc.ncbi.nlm.nih.gov/articles/PMC2654802/)
- [Google Colab documentation](https://research.google.com/colaboratory/faq.html)

See the differential-expression module's source manifests and the reference run's `sessionInfo.txt` for provenance and software versions. Updates will document both new results and remaining methodological limitations.
