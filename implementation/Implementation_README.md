# PTCDiscovery — Implementation Guide

**Version 1.0 · Literature retrieval, LLM integration and local interface**

This module connects the project's differential-expression results to scientific abstracts and presents them through an interactive application. Select a cancer, choose a gene set and ask a question. The application supplies the selected evidence to your chosen language model and returns a sourced explanation.

This is **retrieval-augmented generation (RAG)**: relevant information is supplied at inference time. The language model is not trained or fine-tuned on this dataset.

> A populated `rag_data` directory is required to launch the application. Downloading the prepared data from the repository's **Releases** section is the recommended setup. Rebuilding it from the limma outputs is also possible, but takes longer.

## 1. Repository organization

The source files are grouped by responsibility:

| Location | Contents |
|---|---|
| `Data_Prep/` | `PrepareRAGData.py`, `AnnotateGenes.py`, `CollectLiterature.py` |
| `Answers_Build/` | `LiteratureSearch.py`, `BuildRAGContext.py`, `AskRAG.py`, `LLMProviders.py` |
| `AppRAG.py` | Streamlit application entry point |
| `requirements.txt` | Dependencies for this implementation module |
| `Implementation_README.md` | This guide |
| `rag_data/` | Prepared data downloaded from a Release or generated locally |

Names are case-sensitive on some systems. This guide assumes the parent folder is named `implementation` and the two source folders are named exactly `Data_Prep` and `Answers_Build`.

**Version 1 execution convention:** the current imports and default paths expect the Python scripts to be together. The repository folders organize the source; use the flat execution copy below to run it without changing imports. Merely running `AppRAG.py` from above those folders does not resolve the imports.

## 2. Quick start — recommended

### Install the dependencies

Download or clone the repository, open a terminal at its root, then run:

```bash
cd implementation
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

These commands target macOS/Linux. On Windows PowerShell, activate with `.venv\Scripts\Activate.ps1`. Python 3.12 is a practical starting environment. The requirements file lists direct dependencies rather than a complete version lockfile.

### Create a flat execution directory

From `implementation/`:

```bash
mkdir -p runtime
cp Data_Prep/*.py runtime/
cp Answers_Build/*.py runtime/
cp AppRAG.py runtime/
```

Keep the source folders in the repository. If you later update a script there, copy the updated version into `runtime/` before launching again. Alternatively, manually place all eight Python scripts together in a separate working directory.

### Download the data Release

1. Open the repository's **Releases** section.
2. Download the prepared **RAG Data** archive associated with version 1.
3. Extract it and place the complete `rag_data` folder inside `implementation/runtime/`.

The resulting paths must include:

```text
implementation/runtime/AppRAG.py
implementation/runtime/AskRAG.py
implementation/runtime/BuildRAGContext.py
implementation/runtime/LiteratureSearch.py
implementation/runtime/LLMProviders.py
implementation/runtime/rag_data/retained_genes_annotated.csv
implementation/runtime/rag_data/literature_all_cancers/articles.json
implementation/runtime/rag_data/literature_all_cancers/gene_article_links.csv
implementation/runtime/rag_data/literature_all_cancers/collection_summary.csv
```

The other three preparation scripts are also copied into `runtime/`, but are not executed when the application starts.

**Keep the internal data folders intact. Do not flatten the CSV/JSON files, rename `rag_data`, or leave it nested as `rag_data/rag_data/`.** Only the Python scripts are placed together for execution.

You do not need the original expression matrix, trained classification models or R to use the prepared RAG application.

### Launch

From `implementation/`, with the virtual environment active:

```bash
python -m streamlit run runtime/AppRAG.py --server.address localhost
```

Open **http://localhost:8501** if the browser does not open automatically. Stop the server with **Ctrl+C** in the terminal.

## 3. Using the application

1. Select one of the available cancer comparisons.
2. Choose the genes with the largest absolute expression differences, optionally filtered by direction, or select specific retained gene labels.
3. Choose a documentary topic: biological function, therapeutic research, prognosis/biomarkers or microenvironment.
4. Select a provider, enter the exact model identifier and supply your own API key.
5. Enter a question or use a suggested question, then submit it.
6. Expand the literature, statistical results or usage sections as needed. Download the answer or audit for inspection.

Example: select **LUAD**, choose **FABP4**, and ask: “What do the supplied sources report about FABP4, and how should its lower expression in this analysis be interpreted?”

The selected cancer and genes define the scope. Writing another gene or cancer in the question does not automatically change those controls. The application displays earlier answers for the selected cancer, but **each API request is independent**: previous turns are not included as conversation memory.

The interface supports Gemini, OpenAI, Groq, Anthropic, Mistral, OpenRouter and xAI. A **Compatible OpenAI** option accepts a custom base URL for services implementing Chat Completions. This is not universal support for every model or API protocol. Choose a text-generation model supported by that provider and account. OpenAI models that require a different endpoint are outside this adapter's scope.

For the custom option, use the service's base URL, not its complete `/chat/completions` endpoint. Only enter a trusted endpoint: the API key and selected context are sent there. A local compatible server can use a localhost HTTP address.

### Credentials, data transfer and costs

- No personal API key is embedded in the supplied scripts.
- The key entered in the interface is kept in session memory; it is not intentionally written into answer files, audits or the repository.
- Keys are separated by provider; changing a custom endpoint also changes the associated key field.
- The selected question, statistical context and retrieved abstracts are sent to the chosen provider. Provider-side retention and billing are governed by that service's policies.
- Public data-service URLs are not private credentials. Accessing a public scientific API is different from supplying a personal LLM API key.
- The application does not make an API plan free. Access, quotas and charges depend on the user's account and chosen model.

This version is intended for **local execution**. Publishing the source on GitHub does not host the application. A public multi-user deployment needs its own review of endpoint restrictions, credentials, session isolation and storage.

## 4. What each script does

| Script | Responsibility | Main inputs and outputs |
|---|---|---|
| `PrepareRAGData.py` | Consolidates the per-comparison limma results, checks their structure and separates retained genes from the complete statistical results. | Limma result files → `all_results.csv`, `retained_genes.csv`, `cancer_summary.csv` |
| `AnnotateGenes.py` | Maps original gene labels to HGNC identifiers and approved symbols. Flags unresolved, ambiguous and colliding mappings for review. | Retained genes + HGNC reference → annotated genes, mapping table and review CSV |
| `CollectLiterature.py` | Queries Europe PMC for selected cancer–gene pairs and stores abstracts, metadata, links and collection coverage. | Annotated genes → `literature_all_cancers/` |
| `LiteratureSearch.py` | Searches the locally collected corpus using TF-IDF and cosine similarity after cancer/gene filtering. Excludes flagged publications and specified unsuitable associations. | Local articles + link table → ranked candidate abstracts |
| `BuildRAGContext.py` | Selects statistical results, retrieves relevant abstracts, deduplicates sources and builds a structured evidence context. | Annotated genes + literature search + collection status → context dictionary/JSON |
| `AskRAG.py` | Combines the question, context and scientific instructions, calls the provider adapter, checks citation identifiers and saves the response and audit. | Context + question → answer Markdown and audit JSON |
| `LLMProviders.py` | Translates the common generation request into the provider's REST format and normalizes responses and errors. | Provider/model/key + prompt → generated text, usage and completion status |
| `AppRAG.py` | Provides the Streamlit interface, constructs the request and displays answers, evidence, statistics and downloads. | User selections → calls to the existing retrieval and generation functions |

**Collection and retrieval are different:** `CollectLiterature.py` downloads candidate literature ahead of time. `LiteratureSearch.py` searches that saved collection at question time; it does not perform a new internet search.

`BuildRAGContext.py` does not create the entire knowledge base. It assembles a bounded context from an already prepared statistical table and article collection.

## 5. What is inside `rag_data`?

| Path | Role | Needed to launch? |
|---|---|---|
| `retained_genes_annotated.csv` | Selected differential-expression results with gene identifiers and review flags | **Yes** |
| `literature_all_cancers/articles.json` | Scientific abstracts and bibliographic metadata | **Yes** |
| `literature_all_cancers/gene_article_links.csv` | Candidate cancer–gene–article associations | **Yes** |
| `literature_all_cancers/collection_summary.csv` | Search status for each attempted cancer–gene query | **Yes** |
| `all_results.csv`, `retained_genes.csv`, `cancer_summary.csv` | Statistical preparation outputs and summaries | Not read by the current application; retain for traceability |
| `gene_annotations.csv`, `genes_to_review.csv` | Gene mapping details and review queue | Not read directly by the application |
| `references/` | HGNC reference snapshot and associated provenance | Useful for mapping/reproduction |
| Other collection CSV/JSON files | Selected genes, coverage, collection settings and provenance | Useful for inspection/reproduction |
| `contexts/` | Contexts generated for individual requests | Created automatically; old examples are optional |
| `answers/` | Generated answers and audits | Created automatically |
| `literature_all_cancers/cache/` | Cached collection API responses, if retained | Optional; useful when resuming collection |

The old `literature_LUAD/` pilot collection is not required. The application uses `literature_all_cancers/`.

## 6. Retrieval strategy and evidence boundaries

The statistical ranking comes from the limma table, not from the language model or publication count. Genes are ordered by absolute expression effect, with FDR and original label as tie-breakers. Explicit-gene mode restricts this to the requested genes.

The initial literature collection targeted up to 20 genes per cancer, split between up- and down-regulated genes, with up to five candidate abstracts per query. It is a selected corpus, not an exhaustive review of every retained gene. The collected snapshot reported 664 unique articles and 699 cancer–gene–article links before retrieval-time exclusions. Documentation coverage varies substantially between cancers.

The retriever represents titles and complete abstracts with **TF-IDF word and phrase features**, then uses cosine similarity. It does not currently use neural semantic embeddings or a persistent vector database. The entire abstract is kept as a retrieval unit to preserve experimental context. English topic keywords help search the mainly English corpus alongside the user's question; this is not full multilingual semantic retrieval.

Cancer filtering restricts the association table, but an article can discuss several cancers. An association marked `pending` remains a candidate association, not a verified biological conclusion. A gene may remain in the statistical ranking even if no supporting abstract was collected.

Citation checks verify that cited identifiers were supplied to the model. They **do not automatically prove that each claim faithfully represents the cited article**.

## 7. Rebuild the data instead of downloading it

Use a separate working copy if you want to preserve the released dataset.

1. Obtain the completed limma result folders described in the [differential-expression module](../differential_expression/).
2. Copy the eight implementation scripts together as in the quick start.
3. Set the limma input directory in `PrepareRAGData.py` to your actual result location. Inspect the configuration at the top of each preparation script before executing it; preparation paths are not configured through the interface.
4. Run the preparation stages in order from that working directory:

```bash
python PrepareRAGData.py
python AnnotateGenes.py
python CollectLiterature.py
```

5. Inspect the mapping review file, collection errors and coverage reports.
6. Check that the four mandatory application input files listed above exist.
7. Launch `AppRAG.py` through Streamlit.

Do not rerun preparation blindly over an existing data snapshot. The preparation script may refuse an existing output directory. External databases evolve, so a fresh collection may differ from the released snapshot. Collection time depends on the number of queries, API response time and available cache; no GPU is required for the current retrieval pipeline.

## 8. Errors, usage and reproducibility

| Symptom | Action |
|---|---|
| `ModuleNotFoundError` for a project script | Use the flat execution directory with all eight `.py` files together. |
| Missing CSV/JSON | Check the four required paths and avoid an extra archive-wrapper folder. |
| HTTP 401/403 | Check the provider, key and permissions. |
| HTTP 404 | Check the exact model identifier and endpoint. |
| HTTP 429 | Check request/token limits or account quota. |
| HTTP 503 | Service temporarily unavailable or overloaded; retry later. This does not establish a token-quota violation. |
| No documentation for a gene | Distinguish “not searched”, “no results”, “collection failed” and “filtered out”. |
| Empty or truncated answer | Inspect the completion status, context volume and generation limit. |

The interface exposes each response's reported usage. Provider field names and token accounting differ. The output-token ceiling is not an input-token limit; lowering the number of genes or abstracts reduces the input context. No application-level automatic retry or provider fallback is performed.

Audits retain the question, selected context, scientific instructions, provider/model, returned usage, completion status and citation warnings. Contexts also record input-file hashes. Treat generated answers as unreviewed output until their claims have been checked against the evidence.

## 9. Scientific scope

The application explores associations in bulk gene-expression data and the supplied literature. Differential expression does not establish causation, validate a drug target or justify a treatment. Tissue composition, adjacent-normal controls, cohort imbalance and other confounders affect interpretation. Abstract-only evidence can omit important experimental details.

This module is a research prototype, not a diagnostic or clinical decision system. Future work includes broader literature coverage, semantic or hybrid retrieval, systematic answer evaluation and improved evidence review.
