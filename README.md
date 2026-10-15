# PTCDiscovery

**Version 1.0 — From cancer transcriptomics to literature-grounded exploration**

PTCDiscovery is a bioinformatics research prototype combining gene-expression classification, differential-expression analysis and a retrieval-augmented language-model interface. It connects statistical observations from public cancer data to scientific literature so users can explore genes and potential research directions with inspectable sources.

The project is independently developed and inspired by AI-assisted biomarker and therapeutic-target discovery approaches, including work associated with Sanofi and Owkin. It is not affiliated with, endorsed by or a reproduction of their proprietary platforms.

> **Start here:** to use the interactive application, download the prepared **RAG Data** asset from this repository's **Releases** section and follow the [implementation guide](implementation/Implementation_README.md). The application requires the extracted `rag_data` folder. You do not need to retrain the classifiers or rerun limma.

## What version 1 includes

- Cleaning and inspection of a public TCGA gene-expression dataset, including structured missing-value analysis.
- A logistic-regression classifier for normal versus primary-tumor tissue.
- A PyTorch multilayer perceptron for classification across 32 cancer categories.
- A two-stage classification pipeline evaluated on a held-out patient-grouped test set.
- Anatomical metadata enrichment and cohort preparation for differential expression.
- Separate limma analyses for 11 eligible cancer–tissue comparisons.
- Gene-identifier mapping and collection of candidate scientific abstracts.
- Local literature retrieval and source-linked LLM explanations.
- A Streamlit interface with user-supplied API keys and multiple provider adapters.

The classifier pipeline and RAG application are separate components. In version 1, the RAG interface explores a cancer selected by the user; it does not diagnose a new sample or automatically run the classifiers.

## Run the application

### 1. Get the code and install the implementation environment

Download or clone this repository. From its root:

```bash
cd implementation
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The activation command is for macOS/Linux. On Windows PowerShell, use `.venv\Scripts\Activate.ps1`.

### 2. Put the scripts together for execution

Version 1 organizes the source by responsibility, but the current imports expect the scripts in one working directory:

```bash
mkdir -p runtime
cp Data_Prep/*.py runtime/
cp Answers_Build/*.py runtime/
cp AppRAG.py runtime/
```

This creates a flat execution copy while preserving the repository organization. Refresh these copies after updating the source. You can also manually copy all eight implementation scripts into a single working folder.

### 3. Download RAG Data

Open this repository's **Releases** section and download the version 1 prepared RAG Data archive. Extract it so that this path exists:

```text
implementation/runtime/rag_data/retained_genes_annotated.csv
```

The complete `literature_all_cancers` folder must remain inside `rag_data`. **Keep the data hierarchy intact: place the Python scripts together, but do not move every CSV and JSON into one flat folder.**

The prepared archive is the recommended route. Rebuilding the data from the limma outputs is supported through the preparation scripts and takes longer. Both routes must produce a populated `rag_data` folder before the application can run.

### 4. Launch

From `implementation/`, with the environment active:

```bash
python -m streamlit run runtime/AppRAG.py --server.address localhost
```

Open **http://localhost:8501**. Choose a cancer, select genes, choose an LLM provider and model, enter your API key and ask a question. Results, sources and statistical details can be inspected and downloaded.

The key is held in session memory rather than written to generated files. The selected question and evidence are sent to the chosen provider, whose access conditions, retention policies and charges apply. The application does not provide shared private keys or guarantee free API usage.

For complete installation, folder descriptions, troubleshooting and reconstruction instructions, see [Implementation README](implementation/Implementation_README.md).

## Project modules

| Module | Purpose | Documentation |
|---|---|---|
| Data processing | Expression/phenotype integration, quality checks and feature filtering | Processing scripts and existing module documentation |
| Classification models | Binary logistic regression, multiclass MLP and model persistence | Model scripts and existing module documentation |
| Final evaluation | Two-stage predictions and held-out test metrics | Final-test scripts and existing module documentation |
| [Differential expression](differential_expression/) | Metadata enrichment, cohort selection, limma analysis and result exports | README and analysis summaries in that directory |
| [Implementation](implementation/Implementation_README.md) | Gene annotation, literature collection, retrieval, LLM adapters and interface | Detailed implementation guide |
| Repository Releases | Large prepared datasets and analysis-result archives | Instructions attached to each Release |

Follow the individual module READMEs for upstream input paths, artifacts and execution details. The implementation requirements file covers the RAG module only; upstream model training additionally uses packages such as PyTorch and joblib, while differential-expression analysis uses R and limma.

## Data and analysis strategy

### Public expression data

The expression source is the UCSC Xena TCGA Pan-Cancer Atlas dataset:

[`EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena`](https://xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net)

The supplied expression values are normalized and represented on a `log2(normalized_value + 1)` scale. They are not raw RNA-seq counts. Phenotype annotations provide sample type and cancer-cohort labels. Anatomical enrichment uses information obtained from the [Genomic Data Commons](https://gdc.cancer.gov/).

The working merged dataset contained 11,060 samples. Filtering 219 constant features from 20,531 initial gene features left 20,312. Missing-value profiling identified five sample missingness patterns and 3,588,487 missing entries affecting 3,977 genes. These incomplete genes were removed, retaining **16,335 complete gene features**, without imputing expression values.

Removing all incomplete genes avoided discarding entire heavily affected cohorts, but also removed potentially informative features. The globally performed feature filtering remains a limitation of the predictive evaluation; a stricter future replication should fit preprocessing choices within training data.

### Classification

The classification subset contains **10,438 samples**: 9,701 primary tumors and 737 solid-tissue normal samples. Splits group samples by patient to avoid placing the same patient in training and evaluation sets:

| Partition | Samples |
|---|---:|
| Training | 6,262 |
| Validation | 2,087 |
| Final test | 2,089 |

The first model distinguishes normal from tumor tissue using logistic regression. The second is a feed-forward MLP with layers `16,335 → 128 → 64 → 32`, batch normalization, ReLU and dropout. It classifies tumor samples among the retained cancer categories. It is not a recurrent neural network.

### Differential expression

Samples are selected using cohort membership, sample type and available specimen-site evidence. Primary tumors are compared with compatible solid-tissue normal controls. Insufficiently supported comparisons are excluded rather than generating every possible cancer–organ combination.

The 11 retained comparisons are:

| Comparison | Cancer | Tissue site |
|---|---|---|
| `BRCA_breast` | Breast invasive carcinoma | Breast |
| `COAD_colon` | Colon adenocarcinoma | Colon |
| `KIRC_kidney` | Kidney clear cell carcinoma | Kidney |
| `KIRP_kidney` | Kidney papillary cell carcinoma | Kidney |
| `LIHC_liver` | Liver hepatocellular carcinoma | Liver |
| `LUAD_lung` | Lung adenocarcinoma | Lung |
| `LUSC_lung` | Lung squamous cell carcinoma | Lung |
| `PRAD_prostate` | Prostate adenocarcinoma | Prostate |
| `STAD_stomach` | Stomach adenocarcinoma | Stomach |
| `THCA_thyroid` | Thyroid carcinoma | Thyroid |
| `UCEC_endometrium` | Uterine corpus endometrioid carcinoma | Endometrium |

These are **11 comparisons across nine tissue sites**, not 11 different organs.

The R analysis uses limma on the supplied log-transformed expression values. A tumor-minus-normal contrast estimates expression differences. Patient blocking and `duplicateCorrelation` account for paired observations in eligible comparisons. TREAT tests an effect threshold of 1 on the supplied expression scale, with trend and robust variance moderation; Benjamini–Hochberg FDR is controlled at 0.05 within each comparison.

Effects on `log2(normalized_value + 1)` should not be described as exact raw-expression fold changes. PCA and other plots support quality inspection; PCA does not determine differential-expression significance. See the differential-expression module for parameters, sample counts and output interpretation.

### Literature and RAG

Selected gene labels are mapped to HGNC identifiers. Unresolved or ambiguous mappings are flagged. Europe PMC provides candidate abstracts and publication metadata linked to cancer–gene queries.

At question time, the application selects numerical results from the prepared limma table, retrieves candidate abstracts from the local collection and supplies both to the chosen LLM. The model is asked to distinguish observations, experimental findings and hypotheses and to cite supplied sources.

Version 1 uses **TF-IDF and cosine similarity**, with complete abstracts as retrieval units. It does not use neural embeddings, full-text article retrieval or a continuously updated web search. Numerical gene rankings are computed outside the LLM.

## Reported held-out classification results

The final test set contained 148 normal samples and 1,941 tumors.

| Metric | Result |
|---|---:|
| Binary accuracy | 99.33% |
| Binary balanced accuracy | 98.39% |
| Binary ROC-AUC | 0.9994 |
| Tumors detected and correctly typed by the full chain | 1,863 / 1,941 — 95.98% |
| Tumors detected with the true cancer in the top three | 1,926 / 1,941 — 99.23% |
| Full-chain accuracy, including normal samples | 96.07% |
| Full-chain macro-F1 | 0.9369 |

The binary stage missed 10 tumors and produced four false alerts. The end-to-end tumor metrics include first-stage misses. Performance varies by class; for example, rectum adenocarcinoma had a reported recall of 28.57% in the full-chain test.

These are internal held-out results from the processed TCGA dataset, not evidence of clinical validity or external generalization. LLM provider adapters have been checked with simulated responses; successful use with one provider does not validate every model or every generated claim.

## Limitations and next steps

- Normal controls are generally tissue samples from cancer cohorts, not a population of healthy volunteers.
- Bulk expression reflects cell composition as well as cell-intrinsic regulation. Relevant clinical, technical and biological confounders may remain.
- Differential expression, predictive importance and causality are different concepts. None alone establishes a therapeutic target.
- Cancer classes, control counts and literature coverage are uneven. Statistical absence from the retained table does not mean a gene is unexpressed.
- Abstract-only retrieval and lexical matching can miss evidence or omit experimental context. Candidate literature associations require review.
- LLMs can misinterpret evidence even when citation identifiers are valid. Sources and claims must be checked.
- External validation, calibration, systematic RAG evaluation and broader evidence collection remain future work.

Version 1 is a functional research prototype. Updates will focus on scientific validation, retrieval quality, reproducibility and usability. It is **not intended for diagnosis, treatment selection or clinical decision-making**.

## Data resources

- [UCSC Xena](https://xenabrowser.net/) — expression and phenotype data.
- [NCI Genomic Data Commons](https://gdc.cancer.gov/) — clinical and biospecimen metadata.
- [HGNC](https://www.genenames.org/) — gene nomenclature and identifier mapping.
- [Europe PMC](https://europepmc.org/) — publication records and abstracts.
- [limma](https://bioconductor.org/packages/limma/) — differential-expression modeling.

Retain source identifiers and provenance when reusing the outputs. Scientific publications and source datasets remain subject to their respective access and reuse terms.
