# GTEx target-gene sex comparisons on size-factor-normalized expression.
# Usage: Rscript GTEx_Sex_Welch.R paired|unpaired expression.csv output.csv
# Male minus female; two-sided Welch tests; BH over 20 tests per cohort.
# "Paired" denotes the complete five-region cohort, not a paired sex test.
# Output: one CSV with 20 rows, four target genes by five GTEx regions, reporting
# sample counts, means, medians, the male/female log2 mean ratio, Welch t/df,
# raw p, the 95% CI of the male-minus-female mean difference and the BH20 p.
# Donors contribute several regions and the paired and unpaired cohorts overlap,
# so rows within a cohort are not independent samples of distinct donors.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3 || !args[1] %in% c("paired", "unpaired")) {
  stop("Usage: Rscript GTEx_Sex_Welch.R paired|unpaired expression.csv output.csv")
}
cohort <- args[1]
if (file.exists(args[3])) stop("Output already exists; choose a new output path.")
d <- read.csv(args[2], stringsAsFactors = FALSE, check.names = FALSE)
required <- c("gene", "ensembl_gene_id", "sample_id", "sex", "site", "expression")
if (!all(required %in% names(d))) stop("Missing required long-expression columns.")
genes <- c("SLIT1", "SLIT2", "SLIT3", "ROBO2")
sites <- c("Brain - Frontal Cortex (BA9)", "Brain - Cerebellar Hemisphere",
           "Brain - Hippocampus", "Brain - Substantia nigra", "Brain - Hypothalamus")
if (anyNA(d[required]) || !setequal(d$gene, genes) || !setequal(d$site, sites) ||
    !all(d$sex %in% c("F", "M")) || anyDuplicated(d[c("gene", "sample_id")])) {
  stop("Unexpected targets, regions, missing metadata, sex labels or duplicate sample/gene rows.")
}
d$expression <- as.numeric(d$expression)
if (any(!is.finite(d$expression)) || any(d$expression < 0)) {
  stop("Expected finite, nonnegative normalized expression, not batch-corrected residuals.")
}
expected_n <- if (cohort == "paired") 315L else 962L
if (length(unique(d$sample_id)) != expected_n || nrow(d) != expected_n * 4L ||
    any(table(d$sample_id) != 4L)) stop("Input does not match the selected cohort size.")
sample_meta <- unique(d[c("sample_id", "sex", "site")])
if (nrow(sample_meta) != expected_n) stop("Inconsistent sample metadata across genes.")
sample_meta$donor <- sub("^([-.]?GTEX)[-.]([^- .]+).*", "GTEX-\\2", sample_meta$sample_id)
if (any(!grepl("^GTEX-", sample_meta$donor))) stop("Unrecognized GTEx sample identifiers.")
if (anyDuplicated(sample_meta[c("donor", "site")]) ||
    length(unique(sample_meta$donor)) != (if (cohort == "paired") 63L else 313L) ||
    any(tapply(sample_meta$sex, sample_meta$donor, function(x) length(unique(x))) != 1L)) {
  stop("Unexpected donor roster or inconsistent donor sex.")
}
expected_m <- if (cohort == "paired") rep(45L, 5) else c(153L, 157L, 143L, 101L, 147L)
expected_f <- if (cohort == "paired") rep(18L, 5) else c(56L, 58L, 54L, 38L, 55L)

rows <- list()
for (gene in genes) {
  if (length(unique(d$ensembl_gene_id[d$gene == gene])) != 1L) stop("Ambiguous target identifier.")
  for (i in seq_along(sites)) {
    x <- d[d$gene == gene & d$site == sites[i], ]
    male <- x$expression[x$sex == "M"]
    female <- x$expression[x$sex == "F"]
    if (length(male) != expected_m[i] || length(female) != expected_f[i]) {
      stop("Unexpected sample counts for ", gene, " / ", sites[i])
    }
    test <- t.test(male, female, var.equal = FALSE, alternative = "two.sided")
    rows[[length(rows) + 1L]] <- data.frame(
      gene = gene, ensembl_gene_id = unique(x$ensembl_gene_id), site = sites[i],
      n_male = length(male), n_female = length(female),
      mean_male = mean(male), mean_female = mean(female),
      median_male = median(male), median_female = median(female),
      log2_fc_male_vs_female = log2((mean(male) + 1e-9) / (mean(female) + 1e-9)),
      welch_t = unname(test$statistic), welch_df = unname(test$parameter), p_value = test$p.value,
      ci95_low_mean_diff_male_minus_female = test$conf.int[1],
      ci95_high_mean_diff_male_minus_female = test$conf.int[2])
  }
}
out <- do.call(rbind, rows)
stopifnot(nrow(out) == 20L)
out$p_adj_bh_20_tests <- p.adjust(out$p_value, method = "BH", n = 20L)
out$significant_bh_0_05 <- out$p_adj_bh_20_tests < 0.05
out$cohort <- cohort
dir.create(dirname(args[3]), recursive = TRUE, showWarnings = FALSE)
write.csv(out, args[3], row.names = FALSE, na = "NA")
cat(cohort, ": 20 comparisons;", sum(out$significant_bh_0_05), "BH-significant.\n")
