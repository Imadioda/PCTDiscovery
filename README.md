# PCTDiscovery

**Version 0 | Active development | Independent bioinformatics research project**

A machine learning pipeline using TCGA gene expression data to distinguish normal from tumor tissue and classify 32 cancer types, with planned extensions for biomarker analysis and literature-supported therapeutic target exploration.

## Overview

I am developing this project as a first-year bioinformatics master's student to build experience across the complete research workflow: public data acquisition, quality control, preprocessing, predictive modeling, biological interpretation, and scientific literature retrieval.

Version 0 implements a two-stage prediction pipeline:

1. A logistic regression model distinguishes normal tissue from primary tumor samples.
2. A PyTorch neural network ranks cancer types for samples predicted to be tumors.

The first integrated pipeline has now been evaluated on the reserved internal test set. Biological interpretation, external validation, and literature integration remain future stages.

The inputs are **bulk RNA-seq tissue samples**, not individual cells. A sample classified as normal does not establish that its donor is cancer-free. This is an experimental research prototype, not a clinically validated diagnostic tool.

The project will receive frequent updates as development progresses. Results, implementation choices, and remaining limitations will be documented throughout.

## Inspiration: Sanofi and Owkin

This project is inspired by the Sanofi–Owkin collaboration announced in 2021, which applies AI to oncology research, including biomarker identification, therapeutic target discovery, and treatment-response prediction.

My goal is to explore a related research workflow at an individual project scale, using public cancer transcriptomic data and independently implemented models.

**This project is not affiliated with or endorsed by Sanofi or Owkin.** It does not reproduce their proprietary platform, federated learning infrastructure, or research results. Their collaboration already focused on cancer; this project draws inspiration from its research objectives rather than claiming to transfer a non-oncology system to oncology.

Reference: [Sanofi's official announcement, November 18, 2021](https://www.sanofi.com/en/media-room/press-releases/2021/2021-11-18-06-30-00-2336966).

## Data sources and dimensions

The project uses **The Cancer Genome Atlas (TCGA)** data distributed through **UCSC Xena**, using the TCGA Pan-Cancer Atlas expression collection.

| Input | Local filename |
|---|---|
| Gene expression | `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz` |
| Phenotype annotations | `TCGA_phenotype_denseDataOnlyDownload.tsv` |

The retained annotation fields are `sample`, `sample_type`, and `_primary_disease`. The latter identifies the diagnostic cohort; it can also be present on normal tissue samples from that cohort.

Expression values are provided on the **`log2(norm_value + 1)`** scale after upstream normalization. No additional logarithm is applied. Model-specific standardization is subsequently fitted on training data only.

Sources:

- [UCSC Xena expression dataset](https://xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net).
- [Pan-Cancer Atlas data hub](https://xenabrowser.net/datapages/?host=https%3A%2F%2Fpancanatlas.xenahubs.net): associated datasets and metadata.
- [UCSC Xena team explanation of the expression transformation](https://groups.google.com/g/ucsc-cancer-genomics-browser/c/LWgst7tuo1w).

The source expression matrix has genes in rows and samples in columns. After transposition and merging, **rows represent samples and expression columns represent genes**.

| Processing stage | Samples | Expression features |
|---|---:|---:|
| After merging and duplicate-name handling | 11,060 | 20,531 |
| After exploratory removal of 219 zero-variance features | 11,060 | 20,312 |
| Complete expression matrix used in version 0 | 11,060 | 16,335 |
| Normal/primary-tumor modeling subset | 10,438 | 16,335 |
| Primary-tumor subset for cancer classification | 9,701 | 16,335 |

Including three identifier/annotation columns, the prepared full table has **11,060 rows and 16,338 columns**.

Here, a “gene” refers to an expression feature. Duplicate-name processing retained two profiles as `SLC35E2_1` and `SLC35E2_2`; strict verification and provenance of the duplicate handling remain part of the reproducibility roadmap.

## Missing-data analysis and preprocessing decisions

The initial diagnostic reported approximately **3.97 million missing values**. These were already present before merging. A missing value is not equivalent to zero expression.

Removing 219 zero-variance features during exploration reduced the count to **3,588,487 NaNs**. Examining the exact missing-gene positions revealed only five sample profiles:

| Missingness profile | Samples | Missing genes per sample |
|---|---:|---:|
| 0 | 9,310 | 0 |
| 1 | 33 | 1,224 |
| 2 | 908 | 1,236 |
| 3 | 211 | 3,547 |
| 4 | 598 | 2,805 |

Within each profile, samples have missing values at exactly the same gene positions. The profiles are strongly associated with diagnostic cohorts. Dropping every incomplete sample would remove the entire acute myeloid leukemia cohort and almost all ovarian, stomach, and esophageal samples.

Comparing missing-gene sets across profiles identified:

- **3,977 incomplete features**, representing 19.58% of the 20,312 remaining features.
- **16,335 features observed in every sample**.
- Only five genes missing exclusively in the 33-sample prostate profile: removing those samples would recover just five additional complete features.

### Version 0 decision

The first prediction pipeline uses the **16,335 complete features without imputation**, preserving all 11,060 samples at the cleaning stage. Original data remain available for later analyses of genes excluded from the common feature set. Each model's sample-type restrictions are applied afterward.

This is a documented trade-off, not evidence that excluded genes are biologically unimportant. The exact technical origin of the missingness blocks remains unresolved.

## Models

### Stage 1: normal versus tumor tissue

The binary model uses **9,701 `Primary Tumor` and 737 `Solid Tissue Normal` samples**. Other sample types are outside the current pipeline's scope. This scope follows the annotations rather than a strict solid-cancer definition: lymphomas annotated as `Primary Tumor` are included, whereas blood-derived leukemia samples are excluded.

- Training-only zero-variance filtering and standardization within a scikit-learn pipeline.
- Regularized logistic regression with L-BFGS and `class_weight="balanced"`.
- Validation comparison of `C = 0.01, 0.1, 1.0`; `C = 1.0` selected.
- A decision threshold of 0.5 fixed before the internal test evaluation.
- Pipeline, input gene order, and metadata saved in `.joblib` format.

### Stage 2: cancer type classification

The neural network is trained on **all true primary tumors in the training partition**, irrespective of the binary model's predictions. The target is `_primary_disease`, encoded into **32 classes**.

Architecture: **16,335 → 128 → 64 → 32**.

- Hidden blocks: linear layer, BatchNorm, ReLU, and dropout (0.3 and 0.2).
- A multilayer perceptron (MLP), not a recurrent neural network.
- A separate standardization fitted on the network's training samples only.
- Unweighted `CrossEntropyLoss` applied to logits; softmax used during inference.
- AdamW with `lr=0.001`, `weight_decay=0.0001`, and batch size 64.
- Up to 100 epochs, with early stopping after ten epochs without validation-loss improvement.
- Best weights, BatchNorm statistics, preprocessing parameters, gene order, and class mapping saved in `.pth` format.

### Integrated inference

Samples first pass through the binary model. Those predicted as tumors are passed to the cancer classifier, which returns its top three classes and their softmax scores. Cancer predictions remain empty for samples predicted as normal.

Each model applies its own saved preprocessing. Cancer scores come from the second model alone; they are not combined with the binary score. Calibration has not yet been assessed, and a high score should not be interpreted as established clinical certainty.

## Evaluation protocol

Both models share the same patient-grouped partitions, with approximately 60% training, 20% validation, and 20% testing. Samples from one patient remain in a single partition.

The initial split attempts to preserve normal/tumor proportions, not the distribution of all 32 cancers. Class coverage and sample counts are checked afterward.

| Partition | Normal | Tumor | Total |
|---|---:|---:|---:|
| Training | 442 | 5,820 | 6,262 |
| Validation | 147 | 1,940 | 2,087 |
| Internal test | 148 | 1,941 | 2,089 |

### Validation results

| Model | Metric | Result |
|---|---|---:|
| Binary | Balanced accuracy | 98.30% |
| Binary | ROC-AUC | 0.9876 |
| Binary | Tumors correctly detected | 1,927 / 1,940 |
| Binary | Normal samples correctly recognized | 143 / 147 |
| Cancer classifier alone | Accuracy | 95.88% |
| Cancer classifier alone | Balanced accuracy | 92.75% |
| Cancer classifier alone | Macro-F1 | 92.73% |
| Cancer classifier alone | Top-3 accuracy | 99.64% |

The best neural-network checkpoint was epoch 9; training stopped at epoch 19. Validation performance varied across cancers, with rectal adenocarcinoma showing the lowest recall (23.3%). These validation results were used during model selection.

### Final internal test results — version 0

The integrated pipeline was evaluated on **all 2,089 reserved test samples**, without retraining or adjusting the binary threshold during testing.

| Component | Metric | Result |
|---|---|---:|
| Binary model | Accuracy | 99.33% |
| Binary model | Balanced accuracy | 98.39% |
| Binary model | ROC-AUC | 0.9994 |
| Binary model | Tumors missed | 10 / 1,941 |
| Binary model | False alarms on normal tissue | 4 / 148 |
| Full pipeline, true tumors only | Detected and correctly classified, top-1 | **1,863 / 1,941 (95.98%)** |
| Full pipeline, true tumors only | Detected and correct cancer in top-3 | **1,926 / 1,941 (99.23%)** |
| Full pipeline, all samples | Accuracy: normal or specific cancer | **96.07%** |
| Full pipeline, all samples | Macro-F1 across 33 final labels | **93.69%** |

Binary confusion matrix:

| Actual / Predicted | Normal | Tumor |
|---|---:|---:|
| Normal | 144 | 4 |
| Tumor | 10 | 1,931 |

A tumor blocked by the first model counts as an error in both full-pipeline tumor metrics. Overall final accuracy also includes normal samples correctly classified as normal: **2,007 / 2,089**.

Rectal adenocarcinoma remains the main weakness, with **8 / 28 samples correctly classified (28.57% recall)**. Cholangiocarcinoma has 100% recall on only three samples, but 37.5% precision, illustrating why class-level metrics and sample sizes matter. Random examples are displayed for illustration only; reported metrics use the entire test set.

The test script generates `predictions_test.csv`, `confusion_chaine_complete.csv`, and `rapport_test_final.txt` under `resultats_test_final/`. The latest simplified script evaluates the binary stage and the full chain; it does not separately report standalone network test performance.

### Methodological limitations

- **Global feature filtering preceded the split.** Selection of complete genes and the initial variance filtering used the full dataset. Model fitting used training samples only, but the full preprocessing workflow was not isolated from held-out data; these results remain exploratory.
- A stricter evaluation must fit data-dependent feature selection on training data and specify how missing values in future samples will be handled.
- Upstream public-data preprocessing and cohort, tissue, or batch effects may influence performance. Patient grouping does not remove all such biases.
- Predictions may reflect tissue of origin as well as cancer-specific expression. Some cohorts lack comparable normal tissue controls.
- This is an internal evaluation, not external validation. The test has now been inspected and must not be reused for tuning while still being described as untouched.
- Scores are not yet calibrated, and rare-class estimates are uncertain.
- Predictive importance, biological association, and causal evidence must be distinguished in subsequent interpretation.

## Repository structure

**Proposed GitHub layout — migration is still pending.** Current Python scripts remain at the project root. Moving them into these folders requires updating package imports and data/model paths first.

```text
SanofiLikeCancerDetector/
├── README.md
├── requirements.txt
├── .gitignore
├── preprocessing/
│   ├── __init__.py
│   ├── SampleDataClean.py
│   ├── DataMerge.py
│   └── ExpressionDataClean.py
├── models/
│   ├── __init__.py
│   ├── ModeleBinaire.py       # proposed binary-training filename
│   ├── NeuralNtwr.py
│   ├── NrnModelTrain.py
│   └── PipelineModel.py
├── LancerPrediction.py       # final evaluation entry point
├── data/                     # proposed location; excluded from Git
│   ├── raw/
│   └── processed/
├── resultats_binaire/        # local generated artifacts
│   ├── partitions.csv
│   └── modele_binaire.joblib
├── resultats_cancer/
│   └── modele_cancer.pth
└── resultats_test_final/
```

Existing filenames are kept here to match the implementation. Earlier local copies of the final evaluation script may be named `PipelineTest.py`. Run the evaluation as a regular Python script, not through pytest.

## Installation and execution

Development environment: Python 3.13 in Miniconda on macOS. Exact package versions from the training environment have not yet been pinned; `requirements.txt` lists direct dependencies rather than providing a fully locked environment.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Download the source datasets from UCSC Xena. **For the current flat script layout**, place both input files at the project root using the filenames listed above. The `data/raw/` location is planned for the future refactor.

Execution order:

1. Prepare `dtf_modele` using the preprocessing scripts.
2. Train the binary model and save both `partitions.csv` and `modele_binaire.joblib`.
3. Train the neural network with `NrnModelTrain.py`, reusing the existing partitions.
4. Run `LancerPrediction.py` after fixing the models and decision threshold.

```bash
# From the project root, with source data and binary partitions available:
python NrnModelTrain.py

# With both trained model files available:
python LancerPrediction.py
```

These commands apply to the current flat layout, not the proposed package structure. `NeuralNtwr.py` should only define the network and import its dependencies; it should not load expression data or launch training.

Large matrices and model checkpoints are excluded from Git by default. Unless checkpoints are distributed separately, users must train the models locally before running inference. Reloading checkpoints also requires the matching network architecture and compatible library versions.

## Runtime and resource considerations

Reading the expression file, transposing matrices, copying DataFrames, and analyzing missingness can take substantial time and memory. A dense 11,060 × 20,531 float64 matrix alone occupies approximately **1.82 GB**, excluding pandas overhead and additional copies.

Preprocessing and training may take several minutes or longer depending on hardware. In the current implementation, some imports rerun data preparation, so loading may account for much of the delay before prediction. Caching prepared data is planned.

The evaluation script reports start/end times and total duration. macOS notifications use `osascript` and depend on system permissions. Optional `tqdm` progress bars can monitor loops, but cannot automatically track the internal progress of a single monolithic file-loading operation.

## Progress and roadmap

### Completed in version 0

- [x] Expression/phenotype merging by sample identifier.
- [x] Initial duplicate-name processing and variance diagnostics.
- [x] Sample-level, gene-level, and exact-profile missingness analysis.
- [x] Assessment of cohort losses and overlapping missing-gene sets.
- [x] Construction of a 16,335-feature matrix without imputation.
- [x] Shared patient-grouped partitions.
- [x] Binary logistic regression training, validation, and saving.
- [x] Multiclass PyTorch MLP training, validation, and saving.
- [x] Integration of the two models into an inference pipeline.
- [x] Full-pipeline internal test evaluation and summary reporting.

### Next development stages

- [ ] Refactor the repository into packages, update paths, and cache prepared data.
- [ ] Pin dependency versions and record data provenance and preprocessing details.
- [ ] Verify duplicate handling for `SLC35E2` and investigate upstream missingness and transformations.
- [ ] Rebuild evaluation with training-only feature selection and an appropriate independent assessment.
- [ ] Compare the network with a simpler multiclass baseline and investigate cohort-specific errors.
- [ ] Assess calibration, uncertainty thresholds, and external generalization.
- [ ] Interpret predictive gene signatures and assess their stability.
- [ ] Perform tumor/normal differential expression using a method compatible with normalized expression; do not treat this matrix as raw counts for DESeq2/edgeR.
- [ ] Explore functional enrichment, biological pathways, and molecular subgroups.
- [ ] Prioritize candidate genes using statistical evidence, biological context, and interaction networks.
- [ ] Build a retrieval-augmented generation (RAG) layer connecting genes, proteins, pathways, and documented interventions to scientific literature.
- [ ] Generate source-backed reports from computed results and retrieved evidence.
- [ ] Develop an interface for exploring predictions and biological findings.

The long-term objective is a research prototype that connects predictive modeling to evidence-supported biological hypotheses. Therapeutic candidates would remain research hypotheses requiring further validation, not treatment recommendations.

## Ongoing updates

Version 0 completes the first working predictive pipeline. Biological interpretation, methodological improvements, and literature integration are still in progress. Future updates will document both improvements and remaining limitations.

Raw data, checkpoints, and local output folders are excluded from Git by default. Selected lightweight reports can be published separately, while the README retains the key metrics and their evaluation context.
