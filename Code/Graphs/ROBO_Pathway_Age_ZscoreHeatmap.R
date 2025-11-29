###############################################################################
# UNIVERSAL HEATMAP v7 — **robust, step-by-step commented**
# ---------------------------------------------------------------------------
# MAIN DIFFERENCES (v6 → v7)
#   • Fixed **row** order already existed (row_group_order).
#   • Now also accepts fixed **column** order:
#         col_group_order <- c("20-59", "60-79")
#     — Applies to the label stored in `col_sites`.
#   • Any group not listed is pushed to the end in alphabetical order.
###############################################################################

## ─────────────────────────────── 1. PARAMETERS ───────────────────────────────
file_path <- "HCP ROBO2 Signaling Significant Genes v6.xlsx"    # ← CHANGE HERE

# FIXED ORDERS
row_group_order <- c("Ena/VASP", "LIMK")
col_group_order <- c("20-59", "60-79")       # ← nova linha

# Manual colors (optional)
group_colors_manual <- c("Ena/VASP" = "#4289c7", "LIMK" = "#489c27")
site_colors_manual  <- c("20-59" = "#6accb8", "60-79" = "#f08d91")

# “Anti-artifact” adjustments
eps_sd      <- 1e-4
winsor_pctl <- 0.01

# Column clustering
cluster_cols_global <- FALSE   # FALSE = cluster por grupo; TRUE = global

# Outliers (columns)
remove_outlier_cols <- TRUE
outlier_n_mad       <- 1
# ---------------------------------------------------------------------------

## ─────────────────────────── 2. PACKAGES ─────────────────────────────────────
pkgs <- c("readxl","pheatmap","colorspace","RColorBrewer",
          "openxlsx","matrixStats","dplyr")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install))
  install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

## ─────────────────────────── 3. FUNCTIONS ─────────────────────────────────────
scale_rows_robust <- function(mat, eps = 1e-4, cut = 0.01) {
  lower <- rowQuantiles(mat, probs =  cut,    na.rm = TRUE)
  upper <- rowQuantiles(mat, probs = 1-cut,   na.rm = TRUE)
  mat_w <- pmin(pmax(mat, lower), upper)
  mu    <- rowMeans(mat_w, na.rm = TRUE)
  sd    <- rowSds(mat_w,   na.rm = TRUE)
  sd[sd < eps] <- 1
  out <- sweep(mat_w, 1, mu, "-") |> sweep(1, sd, "/")
  out[!is.finite(out)] <- 0
  out
}
auto_colors <- function(n, type=c("group","site")) {
  type <- match.arg(type)
  if (n <= 12) qualitative_hcl(n, palette = if (type=="group") "Dark3" else "Set3")
  else         rainbow_hcl    (n, c=75, l = if (type=="group") 60      else 75)
}
join_colors <- function(manual, items, type)
  c(manual,
    setNames(auto_colors(length(setdiff(items,names(manual))), type),
             setdiff(items,names(manual))))

## ───────────────────────── 4. READING & CLEANING ────────────────────────────
dados <- read_xlsx(file_path, na = c(""," ","NA","-","–"))
col_1 <- tolower(names(dados)[1])

num_cols <- sapply(dados, is.numeric)
if (col_1 == "group") num_cols[1:2] <- FALSE else num_cols[1] <- FALSE

expr_raw <- as.matrix(dados[, num_cols]); expr_raw[is.na(expr_raw)] <- 0
expr <- if (all(expr_raw >= 0)) log2(expr_raw + 1) else expr_raw

## ─────── 5. ROW GROUPS (Group vs Gene) ────────────────────────────────
if (col_1 == "group") {
  grupos <- dados[[1]]; genes <- dados[[2]]
} else if (col_1 == "gene") {
  genes <- dados[[1]]
  expr_scaled <- scale_rows_robust(expr, eps_sd, winsor_pctl)
  k <- max(2, min(15, round(sqrt(length(genes)/2))))
  clusters <- cutree(hclust(dist(expr_scaled)), k)
  grupos   <- paste0("Cluster_", clusters)
  write.xlsx(data.frame(Gene=genes, Cluster=grupos),
             "gene_clusters.xlsx", overwrite=TRUE)
  message("✔ 'gene_clusters.xlsx' saved (", k, " clusters).")
} else stop("First column must be named 'Group' or 'Gene'.")
rownames(expr) <- make.unique(genes, sep="_")

## ───────────────────── 6. ROBUST Z-SCORE ────────────────────────────────
expr_z <- scale_rows_robust(expr, eps_sd, winsor_pctl)

## ──────── 7. DETECTION & REMOVAL OF OUTLIERS (COLUMNS) ───────────────────
if (remove_outlier_cols && ncol(expr_z) > 2) {
  mean_cor <- rowMeans(cor(expr_z, use="pairwise.complete.obs"))
  cutoff   <- median(mean_cor) - outlier_n_mad * mad(mean_cor)
  out_cols <- names(mean_cor)[mean_cor < cutoff]
  if (length(out_cols)) {
    message("⚠️ Removing ", length(out_cols)," columns: ",
            paste(out_cols, collapse=", "))
    expr_z <- expr_z[ , !(colnames(expr_z) %in% out_cols), drop = FALSE]
  } else message("✔ No column detected as an outlier.")
}

## ───────────────────────── 8. ANNOTATIONS ──────────────────────────────────
col_sites <- gsub("\\.{3}\\d+$", "", colnames(expr_z)) |> trimws()
row_ann   <- data.frame(Group = grupos, row.names = rownames(expr_z))
col_ann   <- data.frame(Site  = col_sites, row.names = colnames(expr_z))

## ─────────────────────── 9. COLORS ────────────────────────────────────────
group_colors <- join_colors(group_colors_manual, unique(grupos),    "group")
site_colors  <- join_colors(site_colors_manual,  unique(col_sites), "site")
ann_colors   <- list(Group = group_colors, Site = site_colors)

## ─────────── 10. ORDER / CLUSTER ROWS ───────────────────────────
desired_rows <- row_group_order
present_rows <- desired_rows[desired_rows %in% unique(grupos)]
extras_rows  <- setdiff(sort(unique(grupos)), desired_rows)
full_row_seq <- c(present_rows, extras_rows)
row_blocks <- lapply(full_row_seq, function(g) {
  matg <- expr_z[grupos == g, , drop = FALSE]
  if (nrow(matg)==0) return(NULL)
  matg[hclust(dist(matg))$order, , drop = FALSE]
})
row_blocks <- row_blocks[!sapply(row_blocks, is.null)]
final_mat  <- do.call(rbind, row_blocks)
row_ann    <- row_ann[rownames(final_mat), , drop = FALSE]

## ─────────── 10b. ORDER / CLUSTER COLUMNS ─────────────────────────
if (!cluster_cols_global) {
  desired_cols <- col_group_order
  present_cols <- desired_cols[desired_cols %in% unique(col_sites)]
  extras_cols  <- setdiff(sort(unique(col_sites)), desired_cols)
  full_col_seq <- c(present_cols, extras_cols)
  col_blocks <- lapply(full_col_seq, function(s) {
    cols <- which(col_sites == s)
    if (length(cols)==0) return(NULL)
    if (length(cols)==1) return(cols)
    submat <- t(expr_z[, cols, drop = FALSE])
    cols[hclust(dist(submat))$order]
  })
  col_blocks <- col_blocks[!sapply(col_blocks, is.null)]
  col_order  <- unlist(col_blocks)
  final_mat  <- final_mat[ , col_order, drop=FALSE]
  col_ann    <- col_ann[col_order, , drop=FALSE]
}

## ───────────────────── 11. FINAL HEATMAP ───────────────────────────────
breaks <- seq(-3, 3, length.out = 256)

pheatmap(final_mat,
         color             = colorRampPalette(rev(brewer.pal(11,"RdBu")))(255),
         breaks            = breaks,
         annotation_row    = row_ann,
         annotation_col    = col_ann,
         annotation_colors = ann_colors,
         cluster_rows      = FALSE,
         cluster_cols      = cluster_cols_global,
         treeheight_row    = 0,
         show_rownames     = TRUE,
         show_colnames     = FALSE,
         fontsize_row      = 14,
         fontsize          = 10,
         main              = "Heatmap",
         filename          = "heatmap_expressao.png",
         border_color      = NA,
         width             = 15,
         height            = 3)

# ──────────────────────────────── END ───────────────────────────────────
# Heatmap saved as "heatmap_expressao.png" (and "gene_clusters.xlsx", if created).