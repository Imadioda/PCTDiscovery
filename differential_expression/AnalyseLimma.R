# AnalyseLimma.R
# Input: outputs of PrepareDEComparisons.py, already on the Xena log scale.
# No raw-count method, voom, extra log transform or StandardScaler is used.
# TREAT tests an effect exceeding MIN_EFFECT on the INPUT log-expression scale.
# Results are exploratory until expression QC and confounding are reviewed.

# ---- Editable configuration --------------------------------------------------
PROJECT_DIR <- "/Users/mac/Desktop/PrepPython2/SanofiLikeCancerDetector"
INPUT_DIR <- file.path(PROJECT_DIR, "resultats_DE_preparation")
OUTPUT_BASE <- file.path(PROJECT_DIR, "resultats_limma")
FDR_THRESHOLD <- 0.05
MIN_EFFECT <- 1.0
MIN_PATIENTS_PER_GROUP <- 10L
MIN_PAIRS_FOR_BLOCK_MODEL <- 3L
USE_TREND <- TRUE
USE_ROBUST <- TRUE

# Optional terminal use: Rscript AnalyseLimma.R /path/to/input /path/to/output
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1L) INPUT_DIR <- cli_args[1L]
if (length(cli_args) >= 2L) OUTPUT_BASE <- cli_args[2L]

# ---- Dependencies ------------------------------------------------------------
needed <- c("limma", "statmod")
missing_packages <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       ". Install them using the commands in the accompanying guide.")
}
if (!is.finite(FDR_THRESHOLD) || FDR_THRESHOLD <= 0 || FDR_THRESHOLD >= 1 ||
    !is.finite(MIN_EFFECT) || MIN_EFFECT < 0) stop("Invalid statistical thresholds.")

write_table <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

read_metadata <- function(path) {
  # Explicit character columns preserve patient IDs and labels.
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
}

save_pdf <- function(path, draw, width = 9, height = 7) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
}

# Use precisely the paths recorded by the Python preparation manifest.
resolve_input <- function(relative_path) {
  if (is.na(relative_path) || !nzchar(relative_path) || grepl("^(/|[A-Za-z]:)", relative_path) ||
      ".." %in% strsplit(gsub("\\\\", "/", relative_path), "/", fixed = TRUE)[[1L]]) {
    stop("Invalid relative path in plan_analyses.csv: ", relative_path)
  }
  full <- file.path(INPUT_DIR, relative_path)
  if (!file.exists(full)) stop("Missing input file: ", full)
  full
}

make_qc <- function(y, samples, out) {
  colors <- ifelse(samples$group == "Normal", "#2474B5", "#C84043")
  qc <- data.frame(
    sample = samples$sample, patient_id = samples$patient_id, group = samples$group,
    mean_log_expression = colMeans(y), median_log_expression = apply(y, 2L, median),
    sd_log_expression = apply(y, 2L, sd), stringsAsFactors = FALSE
  )
  write_table(qc, file.path(out, "qc_echantillons.csv"))
  save_pdf(file.path(out, "qc_distributions.pdf"), function() {
    # No sample is automatically deleted based on these distributions.
    boxplot(y, outline = FALSE, names = FALSE, col = colors, border = colors,
            main = "Sample expression distributions", xlab = "Samples (metadata order)",
            ylab = "Input log2(norm_value + 1)")
    legend("topright", c("Normal", "Tumor"), fill = c("#2474B5", "#C84043"), bty = "n")
  })
  variances <- apply(y, 1L, var)
  selected <- head(order(variances, decreasing = TRUE), min(500L, nrow(y)))
  pca <- prcomp(t(y[selected, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 2L)
  if (ncol(pca$x) < 2L) stop("Insufficient expression rank for a two-dimensional PCA.")
  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  coordinates <- data.frame(sample = samples$sample, patient_id = samples$patient_id,
                            group = samples$group, PC1 = pca$x[, 1L], PC2 = pca$x[, 2L])
  write_table(coordinates, file.path(out, "qc_PCA_coordonnees.csv"))
  save_pdf(file.path(out, "qc_PCA.pdf"), function() {
    plot(coordinates$PC1, coordinates$PC2, col = colors, pch = 19,
         xlab = sprintf("PC1 (%.1f%%)", explained[1L]),
         ylab = sprintf("PC2 (%.1f%%)", explained[2L]),
         main = sprintf("PCA: %d most variable genes (visualization only)", length(selected)))
    legend("topright", c("Normal", "Tumor"), col = c("#2474B5", "#C84043"), pch = 19, bty = "n")
  })
}

analyse_comparison <- function(entry, run_dir) {
  name <- entry$comparison_id
  if (is.na(name) || !grepl("^[A-Za-z0-9_-]+$", name)) stop("Invalid comparison name.")
  out <- file.path(run_dir, name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  samples <- read_metadata(resolve_input(entry$samples_file))
  required <- c("sample", "patient_id", "group", "gdc_project", "_primary_disease")
  if (!all(required %in% names(samples))) stop("Missing required metadata columns.")
  if (anyNA(samples[, required]) || any(samples$sample == "") || any(samples$patient_id == "")) {
    stop("Missing sample, patient or group metadata.")
  }
  if (anyDuplicated(samples$sample)) stop("Duplicate sample identifiers.")
  if (!all(samples$group %in% c("Normal", "Tumor"))) stop("Unexpected group labels.")
  if (anyDuplicated(samples[, c("patient_id", "group")])) {
    stop("Repeated observations for the same patient and condition require review.")
  }
  if (!all(samples$gdc_project == entry$project) ||
      !all(samples[["_primary_disease"]] == entry$disease)) stop("Project/disease mismatch with the plan.")
  table <- read.delim(gzfile(resolve_input(entry$expression_file)),
                      check.names = FALSE, stringsAsFactors = FALSE, row.names = 1L)
  if (!all(vapply(table, is.numeric, logical(1)))) stop("Expression contains non-numeric columns.")
  y <- as.matrix(table)
  rm(table)
  if (anyDuplicated(rownames(y)) || anyDuplicated(colnames(y))) stop("Duplicate expression identifiers.")
  if (!setequal(colnames(y), samples$sample)) stop("Expression and metadata sample sets differ.")
  # Explicit alignment: never assume the two files are ordered identically.
  samples <- samples[match(colnames(y), samples$sample), , drop = FALSE]
  stopifnot(identical(colnames(y), samples$sample))
  if (nrow(y) < 10L || !all(is.finite(y))) stop("Too few genes or nonfinite expression values.")
  if (any(apply(y, 1L, var) == 0)) stop("Constant genes remain; inspect the preparation output.")
  samples$group <- factor(samples$group, levels = c("Normal", "Tumor"))
  n_normal <- sum(samples$group == "Normal")
  n_tumor <- sum(samples$group == "Tumor")
  if (n_normal < MIN_PATIENTS_PER_GROUP || n_tumor < MIN_PATIENTS_PER_GROUP) {
    stop("Fewer than the configured minimum independent patients per group.")
  }
  if (n_normal != as.integer(entry$n_normal) || n_tumor != as.integer(entry$n_tumor) ||
      nrow(y) != as.integer(entry$n_genes_exported)) stop("File dimensions disagree with plan_analyses.csv.")
  normal_ids <- samples$patient_id[samples$group == "Normal"]
  tumor_ids <- samples$patient_id[samples$group == "Tumor"]
  n_pairs <- length(intersect(normal_ids, tumor_ids))
  write_table(samples, file.path(out, "echantillons_utilises.csv"))
  make_qc(y, samples, out)

  # Mean model: one coefficient per condition. Contrast is Tumor minus Normal.
  design <- model.matrix(~ 0 + group, data = samples)
  colnames(design) <- c("Normal", "Tumor")
  rownames(design) <- samples$sample
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient design matrix.")
  write_table(data.frame(sample = rownames(design), design, check.names = FALSE),
              file.path(out, "matrice_design.csv"))
  correlation <- NA_real_
  method <- "independent_patients"
  if (n_pairs > 0L) {
    # Guard against unstable block-correlation estimation from very few pairs.
    # This is a workflow choice, not a universal mathematical threshold.
    if (n_pairs < MIN_PAIRS_FOR_BLOCK_MODEL) {
      stop(sprintf("Only %d paired patients: block-model design needs review. QC files were saved; no independence assumption was substituted.", n_pairs))
    }
    block <- factor(samples$patient_id)
    dc <- limma::duplicateCorrelation(y, design = design, block = block)
    correlation <- dc$consensus.correlation
    if (!is.finite(correlation) || abs(correlation) >= 1) stop("Invalid estimated within-patient correlation.")
    fit <- limma::lmFit(y, design, block = block, correlation = correlation)
    method <- "patient_block_duplicateCorrelation"
  } else {
    fit <- limma::lmFit(y, design)
  }
  contrast <- limma::makeContrasts(Tumor - Normal, levels = design)
  fit <- limma::contrasts.fit(fit, contrast)
  if (any(fit$df.residual <= 0)) stop("Insufficient residual degrees of freedom.")

  # Empirical Bayes TREAT: H0 = |effect| <= MIN_EFFECT.
  # For this preprocessed log2(norm_value + 1) matrix, the effect is a difference
  # on that scale, NOT necessarily log2 of a raw-count expression ratio.
  fit <- limma::treat(fit, lfc = MIN_EFFECT, trend = USE_TREND, robust = USE_ROBUST)
  result <- limma::topTreat(fit, coef = 1L, number = Inf, adjust.method = "BH", sort.by = "none")
  idx <- match(rownames(result), rownames(y))
  result <- data.frame(
    gene = rownames(result), comparison_id = name, cancer = entry$disease,
    project = entry$project, tumor_sites = entry$tumor_sites, normal_sites = entry$normal_sites,
    n_normal = n_normal, n_tumor = n_tumor, n_paired_patients = n_pairs,
    mean_normal_log_expression = rowMeans(y[, samples$group == "Normal", drop = FALSE])[idx],
    mean_tumor_log_expression = rowMeans(y[, samples$group == "Tumor", drop = FALSE])[idx],
    effect_log2_scale = result$logFC, average_log_expression = result$AveExpr,
    moderated_t_treat = result$t, p_value_treat = result$P.Value, fdr_bh_treat = result$adj.P.Val,
    stringsAsFactors = FALSE, row.names = NULL
  )
  if (any(!is.finite(result$p_value_treat)) || any(!is.finite(result$fdr_bh_treat))) {
    stop("Nonfinite test results: inspect model fit.")
  }
  result$selected <- result$fdr_bh_treat <= FDR_THRESHOLD
  result$direction <- ifelse(result$selected,
                              ifelse(result$effect_log2_scale > 0, "Up_in_tumor", "Down_in_tumor"),
                              "Not_selected")
  result$test <- "limma_TREAT"
  result$minimum_effect_tested <- MIN_EFFECT
  result$fdr_threshold <- FDR_THRESHOLD
  result <- result[order(result$fdr_bh_treat, -abs(result$effect_log2_scale)), , drop = FALSE]
  write_table(result, file.path(out, "resultats_complets.csv"))
  write_table(result[result$selected, , drop = FALSE], file.path(out, "genes_retenus.csv"))
  write_table(result[result$direction == "Up_in_tumor", , drop = FALSE], file.path(out, "genes_surexprimes.csv"))
  write_table(result[result$direction == "Down_in_tumor", , drop = FALSE], file.path(out, "genes_sous_exprimes.csv"))

  colors <- ifelse(result$direction == "Up_in_tumor", "#C84043",
                   ifelse(result$direction == "Down_in_tumor", "#2474B5", "#999999"))
  save_pdf(file.path(out, "volcano.pdf"), function() {
    plot(result$effect_log2_scale, -log10(pmax(result$fdr_bh_treat, .Machine$double.xmin)),
         pch = 16, cex = 0.45, col = colors, main = name,
         xlab = "Tumor - Normal effect on input log2 scale", ylab = "-log10(BH-adjusted TREAT p-value)")
    abline(v = c(-MIN_EFFECT, MIN_EFFECT), h = -log10(FDR_THRESHOLD), lty = 2, col = "#555555")
  })
  save_pdf(file.path(out, "MA_plot.pdf"), function() {
    plot(result$average_log_expression, result$effect_log2_scale,
         pch = 16, cex = 0.45, col = colors, main = name,
         xlab = "Mean input log expression", ylab = "Tumor - Normal effect")
    abline(h = c(-MIN_EFFECT, 0, MIN_EFFECT), lty = c(2, 1, 2), col = "#555555")
  })
  save_pdf(file.path(out, "qc_moyenne_variance.pdf"), function() {
    limma::plotSA(fit, main = paste(name, "- mean/variance trend"))
  })
  capture.output(list(
    comparison = as.list(entry), contrast = "Tumor - Normal", method = method,
    n_paired_patients = n_pairs, consensus_correlation = correlation,
    use_trend = USE_TREND, use_robust = USE_ROBUST,
    test = "TREAT: absolute effect exceeds MIN_EFFECT", minimum_effect = MIN_EFFECT,
    fdr_BH_within_comparison = FDR_THRESHOLD,
    scale = "Existing log2(norm_value + 1), no additional transformation",
    limitations = c("No additional age/sex/batch/purity covariates fitted in this version.",
      "Upstream EB++ processing is inherited; downstream modelling cannot undo upstream confounding.",
      "Inspect QC and tissue details before interpreting these exploratory results.",
      "BH correction is within each comparison, not across all cancers.",
      "Normal samples are patient-derived controls, not automatically healthy-donor tissue.",
      "Association does not establish causality or a therapeutic target.")
  ), file = file.path(out, "parametres_et_limites.txt"))
  data.frame(comparison_id = name, status = "completed_exploratory", n_normal = n_normal,
    n_tumor = n_tumor, n_paired_patients = n_pairs, n_genes_tested = nrow(result),
    n_up = sum(result$direction == "Up_in_tumor"), n_down = sum(result$direction == "Down_in_tumor"),
    correlation = correlation, model = method, message = "", stringsAsFactors = FALSE)
}

run_limma <- function() {
  start <- Sys.time()
  plan_path <- file.path(INPUT_DIR, "plan_analyses.csv")
  if (!file.exists(plan_path)) stop("Missing plan: ", plan_path, ". Check INPUT_DIR.")
  plan <- read_metadata(plan_path)
  required <- c("comparison_id", "status", "project", "disease", "tumor_sites", "normal_sites",
                "expression_file", "samples_file", "n_normal", "n_tumor", "n_genes_exported")
  if (!all(required %in% names(plan))) stop("Unexpected preparation manifest schema.")
  if (anyDuplicated(plan$comparison_id)) stop("Duplicate comparison IDs in manifest.")
  work <- plan[plan$status == "exported_for_expression_qc", , drop = FALSE]
  if (!nrow(work)) stop("No exported comparisons found. Inspect Python plan_analyses.csv.")
  # Never overwrite an earlier run or mix old and new results.
  run_dir <- file.path(OUTPUT_BASE, paste0("run_", format(start, "%Y%m%d_%H%M%S")))
  if (dir.exists(run_dir)) stop("Run directory already exists: ", run_dir)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  write_table(plan, file.path(run_dir, "plan_entree.csv"))
  writeLines(c(paste("Start:", start), paste("Input:", normalizePath(INPUT_DIR)),
    paste("FDR:", FDR_THRESHOLD), paste("Minimum effect:", MIN_EFFECT),
    paste("Trend:", USE_TREND), paste("Robust:", USE_ROBUST)), file.path(run_dir, "execution.txt"))
  on.exit(capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt")), add = TRUE)
  summaries <- vector("list", nrow(work))
  for (i in seq_len(nrow(work))) {
    entry <- work[i, , drop = FALSE]
    cat(sprintf("\n[%d/%d] %s\n", i, nrow(work), entry$comparison_id))
    warnings_seen <- character()
    summaries[[i]] <- tryCatch(
      withCallingHandlers(analyse_comparison(entry, run_dir), warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) {
        cat("Not completed:", conditionMessage(e), "\n")
        data.frame(comparison_id = entry$comparison_id, status = "error_or_design_review_required",
          n_normal = NA_integer_, n_tumor = NA_integer_, n_paired_patients = NA_integer_,
          n_genes_tested = NA_integer_, n_up = NA_integer_, n_down = NA_integer_,
          correlation = NA_real_, model = "", message = conditionMessage(e), stringsAsFactors = FALSE)
      }
    )
    if (length(warnings_seen)) {
      writeLines(unique(warnings_seen), file.path(run_dir, paste0(entry$comparison_id, "_warnings.txt")))
      cat("Warnings recorded; review", entry$comparison_id, "_warnings.txt\n")
      summaries[[i]]$message <- paste(c(summaries[[i]]$message, unique(warnings_seen)), collapse = " | ")
    }
    write_table(do.call(rbind, summaries[seq_len(i)]), file.path(run_dir, "synthese_analyses.csv"))
    if (summaries[[i]]$status == "completed_exploratory") {
      cat(sprintf("Completed: %d up, %d down, %d genes tested.\n",
          summaries[[i]]$n_up, summaries[[i]]$n_down, summaries[[i]]$n_genes_tested))
    }
    invisible(gc())
  }
  final <- do.call(rbind, summaries)
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "mins"))
  cat(sprintf("\nFinished in %.1f min. Completed: %d/%d.\nResults: %s\n", elapsed,
              sum(final$status == "completed_exploratory"), nrow(final), run_dir))
  cat("Inspect synthese_analyses.csv, warnings and QC before interpretation.\n")
  invisible(final)
}

run_limma()
