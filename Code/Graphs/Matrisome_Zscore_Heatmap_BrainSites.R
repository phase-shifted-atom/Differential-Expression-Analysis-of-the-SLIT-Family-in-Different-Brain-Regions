# Matrisome Z-score heatmap across 8 brain sites
# -------------------------------------------------------------------
# This script:
#   1. Reads a matrix of Z-scored expression values for Matrisome genes
#      (rows = genes, columns = brain sites).
#   2. Applies an optional winsorization step to reduce the impact of
#      extreme values on the color scale.
#   3. Builds row annotations for Matrisome groups and column annotations
#      for brain sites.
#   4. Plots a  heatmap (no clustering on rows, optional
#      clustering on columns) and saves it as a PNG file.
#
# Expected input file structure ("Limma Matrisome HCP Significant Wilcoxon v2.xlsx"):
#   - Column 1: Matrisome class for each gene (e.g., ECM Glycoproteins, Collagens).
#   - Column 2: Gene symbol.
#   - Columns 3+: Z-score expression values for each brain site.
#
# Note: This script is designed to be fully reproducible. If you change
#       the input file name or the number/order of brain sites, you may
#       need to adjust `site_levels` and the expected column layout.

library(readxl)
library(pheatmap)
library(RColorBrewer)

# ────────────────── 1. Data import and basic preprocessing ──────────────────
dados  <- read_xlsx("Limma Matrisome HCP Significant Wilcoxon v2.xlsx")

groups <- dados[[1]]        # Matrisome classes for each gene
genes  <- dados[[2]]        # gene symbols
expr   <- as.matrix(dados[, 3:ncol(dados)])   # pre-computed Z-scores (one column per brain site)
rownames(expr) <- genes
expr[is.na(expr)] <- 0      # replace any missing Z-scores with 0 to avoid NA issues in pheatmap

## (optional) winsorize extremes to avoid color saturation in the heatmap
## This keeps very large/small Z-scores from dominating the color scale.
winsorize <- function(mat, p = 0.01) {
  qs <- quantile(mat, c(p, 1 - p), na.rm = TRUE)  # lower/upper quantile cutoffs (e.g. 1% and 99%)
  mat[mat < qs[1]] <- qs[1]
  mat[mat > qs[2]] <- qs[2]
  mat
}
expr_w <- winsorize(expr, p = 0.01)  # trim lower/upper 1% of values (robust to outliers)

# ────────────────── 2. Row annotation (Matrisome groups) ──────────────────
# Fix the plotting order of Matrisome groups so that rows are grouped by biology, not by alphabetical order
desired_group_order <- c(
  "ECM Glycoproteins", "Proteoglycans", "Collagens",
  "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors"
)
row_annotation <- data.frame(
  Group = factor(groups, levels = desired_group_order),
  row.names = genes
)
# Reorder rows of the expression matrix to match the desired Matrisome group order
ord <- order(row_annotation$Group)
expr_w          <- expr_w[ord, ]
row_annotation  <- row_annotation[ord, , drop = FALSE]

# ────────────────── 3. Column annotation (brain sites) ──────────────────
col_sites       <- colnames(expr_w)
col_sites_clean <- gsub("\\.{3}\\d+$", "", col_sites)  # remove trailing "...1", "...2" artifacts from duplicated Excel column names
# Define the expected order of brain sites (used as factor levels for consistent column ordering and annotation)
site_levels     <- c("CER","SUB","HTL","FRO","HCP","PAR","TEMP","OCC")
annotation_col  <- data.frame(
  Site = factor(col_sites_clean, levels = site_levels),
  row.names = col_sites
)

# ────────────────── 4. Color palette and value range ──────────────────
limit   <- max(abs(expr_w))          # keep the color scale symmetric around 0 (appropriate for Z-scores)
breaks  <- seq(-limit, limit, length.out = 256)
palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)

# Set distinct colors for Matrisome groups (rows) and brain sites (columns)
# The group colors are chosen to be visually distinct from the site colors to avoid confusion in the legend.
ann_colors <- list(
  Group = c(
    "ECM Glycoproteins"       = "#C187D0",  # light lilac (does not clash with FRO site color)
    "Proteoglycans"           = "#EA82A8",  # light raspberry (distinct from HTL site color)
    "Collagens"               = "#AEA5C5",  # grayish lavender
    "ECM-affiliated Proteins" = "#A1AFB6",  # light slate blue, relatively neutral (far from CER site color)
    "ECM Regulators"          = "#CBA48B",  # light tan (far from HCP/OCC site colors)
    "Secreted Factors"        = "#AF9D96"   # light taupe
  ),
  Site = c(
    CER  = "#83cceb", SUB  = "#b5e6a2", HTL  = "#ff8ad8",
    FRO  = "#7a81ff", HCP  = "#ff9300", PAR  = "#FFD700",
    TEMP = "#76C7C0", OCC  = "#FB8072"
  )
)

# ────────────────── 5. Heatmap plotting ──────────────────
# Plot a publication-ready heatmap:
#   - No clustering on rows (keeps Matrisome groups contiguous).
#   - Clustering on columns (allows unsupervised grouping of brain sites).
#   - Row and column annotations show Matrisome group and brain site, respectively.
pheatmap(
  expr_w,
  scale             = "none",
  color             = palette,
  breaks            = breaks,
  annotation_row    = row_annotation,
  annotation_col    = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows      = FALSE,  # preserve biological grouping of Matrisome classes
  cluster_cols      = TRUE,   # allow hierarchical clustering of brain sites
  treeheight_row    = 0,
  treeheight_col    = 50,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  fontsize          = 10,
  main              = "Matrisome Heatmap Across 8 Brain Sites",
  filename          = "heatmap_expressao_v5.png",  # output PNG file (edit name/path as needed)
  width             = 12,
  height            = 10,
  annotation_names_row = FALSE,
  border_color        = NA
)