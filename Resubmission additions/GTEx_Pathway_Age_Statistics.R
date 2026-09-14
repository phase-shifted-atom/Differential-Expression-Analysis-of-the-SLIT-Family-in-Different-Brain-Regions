# Grouped-age statistics for the 45-gene HCP pathway panel.
# Usage: Rscript GTEx_Pathway_Age_Statistics.R input.xlsx output.csv
# Input: first row contains age labels; first column contains gene symbols.
# The original panel has 197 samples: 114 aged 60-79 and 83 aged 20-59.
# Three tests per gene; BH correction within each 45-gene test family.
# Effects: older minus younger; log2FC is log2(mean older / mean younger).
# Output: one CSV with 135 rows, 45 genes by three tests, reporting raw and
# BH-adjusted p, the location parameter, SE, 95% CI, effect size and log2FC.
# Effect sizes are test-specific: Cohen's d (Student), Glass's delta on the
# younger-group SD (Welch), standardized Wilcoxon r (Mann-Whitney).
# All 197 input samples enter every test; no display filter is applied here.

suppressPackageStartupMessages({ library(readxl); library(effsize); library(rstatix) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript GTEx_Pathway_Age_Statistics.R input.xlsx output.csv")
if (file.exists(args[2])) stop("Output already exists; choose a new output path.")
sheets <- excel_sheets(args[1])
if (length(sheets) != 1L) stop("Expected the original single-sheet pathway input.")
d <- read_excel(args[1], sheet = sheets[1], col_names = FALSE, .name_repair = "minimal")
labels <- trimws(as.character(unlist(d[1, -1], use.names = FALSE)))
genes <- as.character(d[-1, 1][[1]])
if (length(genes) != 45L || anyNA(genes) || anyDuplicated(genes) ||
    length(labels) != 197L || sum(labels == "60-79") != 114L || sum(labels == "20-59") != 83L) {
  stop("Expected 45 unique genes and the original 114 older / 83 younger samples.")
}
matrix <- as.matrix(d[-1, -1])
storage.mode(matrix) <- "double"
if (any(!is.finite(matrix)) || any(matrix < 0)) stop("Expected finite nonnegative pathway expression.")

rows <- list()
for (i in seq_along(genes)) {
  older <- matrix[i, labels == "60-79"]
  younger <- matrix[i, labels == "20-59"]
  student <- t.test(older, younger, var.equal = TRUE)
  welch <- t.test(older, younger, var.equal = FALSE)
  # Preserve the source's NA treatment if a Wilcoxon confidence-interval warning occurs.
  mann <- tryCatch(wilcox.test(older, younger, conf.int = TRUE, exact = FALSE),
                   warning = function(w) list(p.value = NA, estimate = NA, conf.int = c(NA, NA)),
                   error = function(e) list(p.value = NA, estimate = NA, conf.int = c(NA, NA)))
  groups <- data.frame(group = factor(c(rep("60-79", length(older)), rep("20-59", length(younger)))),
                       value = c(older, younger))
  # wilcox_effsize returns the standardized Wilcoxon effect r, not rank-biserial correlation.
  effects <- c(cohen.d(older, younger)$estimate,
               if (sd(younger) == 0) NA else (mean(older) - mean(younger)) / sd(younger),
               wilcox_effsize(groups, value ~ group)$effsize)
  tests <- list(Student = student, Welch = welch, `Mann-Whitney` = mann)
  for (j in seq_along(tests)) {
    test <- tests[[j]]
    difference <- if (j < 3L) unname(test$estimate[1] - test$estimate[2]) else unname(test$estimate)
    rows[[length(rows) + 1L]] <- data.frame(
      Sheet = sheets[1], Gene = genes[i], Test = names(tests)[j], P_Value = test$p.value,
      Adjusted_P_Value = NA_real_, Location_Parameter = difference,
      SE_Difference = if (j < 3L) test$stderr else NA_real_,
      CI_Lower = test$conf.int[1], CI_Upper = test$conf.int[2], Effect_Size = effects[j],
      Log2_FC = if (mean(younger) == 0) NA_real_ else log2(mean(older) / mean(younger)))
  }
}
out <- do.call(rbind, rows)
for (test in unique(out$Test)) {
  selected <- out$Test == test
  stopifnot(sum(selected) == 45L)
  out$Adjusted_P_Value[selected] <- p.adjust(out$P_Value[selected], method = "BH", n = 45L)
}
stopifnot(nrow(out) == 135L)
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
write.csv(out, args[2], row.names = FALSE, na = "NA")
cat("45 genes; 197 input samples; 135 test rows;",
    sum(out$Test == "Welch" & out$Adjusted_P_Value < 0.05, na.rm = TRUE), "Welch BH-significant.\n")
