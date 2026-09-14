# Region-specific GTEx sex sensitivity analysis on frozen integer raw counts.
# Usage: Rscript GTEx_Sex_DESeq2.R paired|unpaired counts.tsv metadata.csv output_dir
# The supplied count matrices retain the original full eligible gene background.
# Each region is fitted independently with ~ sex; the contrast is male vs female.
# All-gene DESeq2 padj and target-panel BH20 are separate statistical families.
# Outputs written to output_dir: deseq2_target_genes_per_site_BH20.csv (20 rows),
# deseq2_per_site_all_genes_combined.csv (all fitted genes for the five regions),
# deseq2_per_site_size_factors.csv and session_info.txt.
# The combined all-gene table is the deposit object for the complete sex models;
# the supplementary workbooks hold only the 20-row target panel per cohort.

suppressPackageStartupMessages({ library(data.table); library(DESeq2) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4 || !args[1] %in% c("paired", "unpaired")) {
  stop("Usage: Rscript GTEx_Sex_DESeq2.R paired|unpaired counts.tsv metadata.csv output_dir")
}
cohort <- args[1]
outdir <- args[4]
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty; choose a new directory.")
}
genes <- c("SLIT1", "SLIT2", "SLIT3", "ROBO2")
sites <- c("Brain - Cerebellar Hemisphere", "Brain - Frontal Cortex (BA9)",
           "Brain - Hippocampus", "Brain - Hypothalamus", "Brain - Substantia nigra")
metadata <- read.csv(args[3], stringsAsFactors = FALSE, check.names = FALSE)
required <- c("sample_id_gct", "sample_id_dots", "donor_id", "sex", "site_label", "SMAFRZE", "SMTS", "AGE")
if (!all(required %in% names(metadata)) || anyNA(metadata[required])) stop("Incomplete sample metadata.")
expected_samples <- if (cohort == "paired") 315L else 962L
expected_donors <- if (cohort == "paired") 63L else 313L
if (nrow(metadata) != expected_samples || length(unique(metadata$donor_id)) != expected_donors ||
    anyDuplicated(metadata$sample_id_gct) || anyDuplicated(metadata$sample_id_dots) ||
    anyDuplicated(metadata[c("donor_id", "site_label")]) || !setequal(metadata$site_label, sites) ||
    !all(metadata$sex %in% c("F", "M")) || !all(metadata$SMAFRZE == "RNASEQ") ||
    !all(metadata$SMTS == "Brain") || !all(metadata$AGE %in% c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79"))) {
  stop("Metadata does not match the original eligible cohort.")
}
if (cohort == "paired" && (!"FIVEPAIRED" %in% names(metadata) ||
    !all(as.character(metadata$FIVEPAIRED) == "TRUE"))) stop("Expected the complete five-region cohort.")
if (any(gsub("-", ".", metadata$sample_id_gct, fixed = TRUE) != metadata$sample_id_dots) ||
    any(tapply(metadata$sex, metadata$donor_id, function(x) length(unique(x))) != 1L) ||
    any(tapply(metadata$AGE, metadata$donor_id, function(x) length(unique(x))) != 1L)) {
  stop("Inconsistent sample identifiers or donor-level metadata.")
}
metadata$sex <- factor(metadata$sex, levels = c("F", "M"))
rownames(metadata) <- metadata$sample_id_dots

raw <- fread(args[2], sep = "\t", check.names = FALSE)
# The paired export uses R-safe sample names; the unpaired export preserves GCT names.
id_col <- if (cohort == "paired") "ensembl_gene_id" else "Name"
symbol_col <- if (cohort == "paired") "gene_name" else "Description"
sample_cols <- if (cohort == "paired") metadata$sample_id_dots else metadata$sample_id_gct
if (!setequal(names(raw), c(id_col, symbol_col, sample_cols))) stop("Raw-count columns do not match the manifest.")
versioned_ids <- as.character(raw[[id_col]])
ids <- if (cohort == "paired") versioned_ids else sub("\\.[0-9]+$", "", versioned_ids)
if (anyNA(ids) || anyDuplicated(ids)) stop("Missing or duplicate Ensembl identifiers.")
annotation <- data.frame(ensembl_gene_id = ids, ensembl_gene_id_gct = versioned_ids,
                         gene_name = as.character(raw[[symbol_col]]))
target_map <- annotation[!is.na(annotation$gene_name) & annotation$gene_name %in% genes, ]
if (nrow(target_map) != 4L || !setequal(target_map$gene_name, genes)) stop("Ambiguous target-gene mapping.")
expected_genes <- if (cohort == "paired") 56200L else 54738L
if (nrow(raw) != expected_genes) stop("Raw counts must retain all ", expected_genes, " genes in the original cohort manifest.")
count_matrix <- as.matrix(raw[, ..sample_cols])
if (!is.numeric(count_matrix) || any(!is.finite(count_matrix)) || any(count_matrix < 0) ||
    any(count_matrix > .Machine$integer.max) || any(abs(count_matrix - round(count_matrix)) > 1e-8)) {
  stop("DESeq2 requires finite nonnegative integer raw counts.")
}
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- ids
colnames(count_matrix) <- metadata$sample_id_dots
rm(raw)
expected_m <- if (cohort == "paired") rep(45L, 5) else c(157L, 153L, 143L, 147L, 101L)
expected_f <- if (cohort == "paired") rep(18L, 5) else c(58L, 56L, 54L, 55L, 38L)
for (i in seq_along(sites)) {
  m <- metadata[metadata$site_label == sites[i], ]
  if (sum(m$sex == "M") != expected_m[i] || sum(m$sex == "F") != expected_f[i]) {
    stop("Unexpected sex counts in ", sites[i])
  }
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tables <- factors <- list()
for (region in sites) {
  cat("Fitting", cohort, "/", region, "\n")
  selected <- metadata$site_label == region
  m <- metadata[selected, , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(count_matrix[, selected, drop = FALSE], m, design = ~ sex)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, quiet = TRUE, parallel = FALSE)
  # Preserve the original default independent-filtering threshold (alpha = 0.1).
  result <- as.data.frame(results(dds, contrast = c("sex", "M", "F"), alpha = 0.1))
  a <- annotation[match(rownames(result), annotation$ensembl_gene_id), ]
  result <- cbind(data.frame(model = "per_site_sex_M_vs_F", site_label = region,
                            coefficient = "sex_M_vs_F", n_male = sum(m$sex == "M"),
                            n_female = sum(m$sex == "F")), a, result)
  result$cohort <- cohort
  tables[[region]] <- result
  factors[[region]] <- data.frame(sample_id = names(sizeFactors(dds)), site_label = region,
                                  sex = as.character(m$sex), size_factor = as.numeric(sizeFactors(dds)))
  rm(dds)
  invisible(gc())
}
all_results <- rbindlist(tables)
targets <- all_results[gene_name %in% genes]
if (nrow(targets) != 20L || anyDuplicated(targets[, .(gene_name, site_label)]) || anyNA(targets$pvalue)) {
  stop("Expected 20 unique estimable target/site tests.")
}
targets$target_bh20_padj <- p.adjust(targets$pvalue, method = "BH", n = 20L)
targets$target_bh20_significant_0_05 <- targets$target_bh20_padj < 0.05
targets <- targets[order(match(gene_name, genes), match(site_label, sites))]
fwrite(targets, file.path(outdir, "deseq2_target_genes_per_site_BH20.csv"), na = "NA")
fwrite(all_results, file.path(outdir, "deseq2_per_site_all_genes_combined.csv"), na = "NA")
fwrite(rbindlist(factors), file.path(outdir, "deseq2_per_site_size_factors.csv"), na = "NA")
writeLines(c(paste("Cohort:", cohort), capture.output(sessionInfo())), file.path(outdir, "session_info.txt"))
cat("20 target comparisons;", sum(targets$target_bh20_significant_0_05), "target-BH significant.\n")
