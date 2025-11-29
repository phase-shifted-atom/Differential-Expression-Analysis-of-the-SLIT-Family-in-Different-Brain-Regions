# ---------------------------------------------------------------------------
# Script: GTEx_Age_Lollipop_ByMatrisomeType.R
# Purpose:
#   Generate a lollipop plot summarising age-related differential expression
#   of hippocampal matrisome genes grouped by Matrisome type.
#
# Input:
#   - "Significant HCP Expression.xlsx"
#       Expected to contain at least:
#         * gene                : gene symbol (e.g. "SLIT1", "COL1A1")
#         * Matrisome_type      : Matrisome class/category
#         * Fishers_combined_p  : Fisher's combined p-value for age effect
#
# Output:
#   - "LollipopChart_Genes_by_Matrisome_Type_v2.png"
#       Publication-ready lollipop plot for a manuscript figure panel, with
#       genes on the y-axis and –log10(Fisher's p) on the x-axis.
#
# Notes:
#   - The script:
#       * Standardises column names (spaces → underscores) for robust piping.
#       * Orders Matrisome types in a biologically meaningful sequence.
#       * Highlights SLIT1/2/3 genes visually, as they are of special interest.
#       * Uses HTML formatting (via ggtext) to italicise gene names on the y-axis.
#   - This script assumes that data pre-filtering (e.g. significance thresholds)
#     has been performed upstream; it focuses on figure generation.
# ---------------------------------------------------------------------------

# Load necessary libraries
# - ggplot2 / ggtext : plotting and rich axis text
# - readxl          : import Excel spreadsheets
# - dplyr           : data wrangling
# - scales          : transformations and helpers (if needed)
library(ggplot2)
library(readxl)
library(dplyr)
library(scales)
library(ggtext)  # for colored axis text

# Read the table of significant HCP genes
# (This file should already contain only the genes to be plotted in the figure.)
data <- read_excel("Significant HCP Expression.xlsx")

# Standardize column names: replace spaces with underscores
# This makes column access safer inside dplyr/ggplot pipelines and avoids
# backticks in code (e.g. "Matrisome type" → "Matrisome_type").
names(data) <- gsub(" ", "_", names(data))

# Create a unique identifier combining gene and Matrisome_type
# This factor will be used on the y-axis to ensure that each gene × class
# combination is treated as a distinct row in the plot.
data <- data %>%
  mutate(gene_cluster = paste(gene, Matrisome_type, sep = "_"))

# Decide which subset of genes to plot
# By default, we use the full set of significant genes imported above.
# If needed for layout or readability, the chunk below can be used to
# restrict to the top N genes by Fisher's combined p-value.
# data_top <- data %>%
#   arrange(desc(Fishers_combined_p)) %>%
#   head(50)
data_top <- data

# Compute –log10(Fisher's combined p-value)
# This transforms very small p-values into larger positive numbers that are
# easier to visualise and compare across genes.
data_top <- data_top %>%
  mutate(negLogP = -log10(Fishers_combined_p))

# 1. Define the exact Matrisome_type order
# This ensures that Matrisome classes appear in a biologically meaningful,
# consistent order in the legend and any grouped summaries.
desired_order <- c(
  "ECM Glycoproteins",
  "Proteoglycans",
  "Collagens",
  "ECM-affiliated Proteins",
  "ECM Regulators",
  "Secreted Factors"
)

# 2. Convert Matrisome_type to an ordered factor
# Using 'desired_order' as levels guarantees that ggplot will respect this
# order, rather than alphabetically sorting the classes.
data_top$Matrisome_type <- factor(
  data_top$Matrisome_type,
  levels = desired_order
)

# 3. Arrange rows by Matrisome_type and then by p-value
# Within each Matrisome class, genes will appear ordered from most to least
# significant (smallest to largest p-value).
data_top <- data_top %>%
  arrange(Matrisome_type, Fishers_combined_p)

# Define the plotting order for gene_cluster
# We reverse the arranged order so that the most significant genes appear
# towards the top of the y-axis (a common convention in lollipop plots).
levels_order <- rev(data_top$gene_cluster)
data_top$gene_cluster <- factor(data_top$gene_cluster, levels = levels_order)

# Determine the left x-axis limit for the baseline segments
# Using the floored minimum ensures that all segments start from a common
# reference value just below the smallest –log10(p) in the dataset.
x_min <- floor(min(data_top$negLogP, na.rm = TRUE))

# Define palettes for Matrisome_type (fills and strokes)
# Colors are chosen to be visually distinct yet harmonious across classes.
# Fills and strokes are separated so that we can:
#   - use softer fills for the point interiors
#   - use more saturated strokes for point outlines and the legend keys.
matrisome_fills <- c(
  "ECM Glycoproteins"       = "#C187D0",  # lilás claro
  "Proteoglycans"           = "#EA82A8",  # framboesa clara
  "Collagens"               = "#AEA5C5",  # lavanda acinzentada
  "ECM-affiliated Proteins" = "#A1AFB6",  # azul-ardósia claro
  "ECM Regulators"          = "#CBA48B",  # tan claro
  "Secreted Factors"        = "#AF9D96"   # taupe claro
)
matrisome_strokes <- c(
  "ECM Glycoproteins"       = "#8E24AA",  # roxo vívido
  "Proteoglycans"           = "#D81B60",  # carmim
  "Collagens"               = "#6B5B95",  # ametista
  "ECM-affiliated Proteins" = "#546E7A",  # azul-cinza médio
  "ECM Regulators"          = "#A05A2C",  # cobre
  "Secreted Factors"        = "#6D4C41"   # mocha
)

# Prepare custom y-axis labels
# We strip the Matrisome_type suffix so that only gene symbols are shown,
# then use HTML/markdown formatting (ggtext) to italicise gene names and
# highlight SLIT1/2/3 in red (genes of particular biological interest).
base_labels <- sub("_.*", "", levels(data_top$gene_cluster))

colored_labels <- ifelse(
  base_labels %in% c("SLIT1", "SLIT2", "SLIT3"),
  paste0("<span style='color:red;'><em>", base_labels, "</em></span>"),
  paste0("<em>", base_labels, "</em>")
)

# Build the lollipop plot
# Each gene (y) is represented by a horizontal segment from x_min to its
# –log10(p) value, with a point colored by Matrisome type at the end of the
# segment. SLIT genes receive an extra red outline for emphasis.
p2 <- ggplot(data_top, aes(
  x = negLogP,
  y = gene_cluster,
  color = Matrisome_type
)) +
  # Grey baseline segments: encode the magnitude of –log10(p) as segment length
  # starting from a common baseline (x_min).
  geom_segment(aes(
    x    = x_min,
    xend = negLogP,
    y    = gene_cluster,
    yend = gene_cluster
  ), color = "gray70") +
  # Main lollipop heads: filled circles colored by Matrisome class
  # (fill = Matrisome_type, color = Matrisome_type).
  geom_point(aes(fill = Matrisome_type, color = Matrisome_type), shape = 21, size = 3, stroke = 0.9) +
  # Additional red circles around SLIT1/2/3
  # This provides a visual cue for genes that are the primary focus of the
  # manuscript, without altering their underlying data values.
  geom_point(
    data = subset(data_top, gene %in% c("SLIT1", "SLIT2", "SLIT3")),
    aes(x = negLogP, y = gene_cluster),
    shape  = 21,
    fill   = NA,
    color  = "red",
    size   = 5,
    stroke = 1
  ) +
  # apply manual fill and color scales
  scale_fill_manual(values = matrisome_fills) +
  scale_color_manual(values = matrisome_strokes) +
  # Constrain x-axis to start at x_min so all baseline segments are visible
  # and the plot does not extend unnecessarily to the left.
  scale_x_continuous(limits = c(x_min, NA)) +
  # Apply HTML-formatted y-axis labels (italic and selectively colored)
  # handled by ggtext::element_markdown in the theme below.
  scale_y_discrete(labels = colored_labels) +
  # axis and legend labels
  labs(
    x    = "-log10(Fisher's Combined p-value)",
    y    = "Gene",
    fill = "Matrisome Type"
  ) +
  # Use a minimal theme and adjust axis/legend text
  # element_markdown enables rendering of the HTML labels created above.
  theme_minimal() +
  theme(
    axis.text.y    = element_markdown(size = 8),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(ncol = 3), color = "none")

# Display the plot in the current device and save a high-resolution PNG file
# for inclusion in the manuscript. Width/height are tuned for legibility of
# gene labels and legend in a vertical layout.
print(p2)
ggsave("LollipopChart_Genes_by_Matrisome_Type_v2.png", plot = p2,
       width = 7.3, height = 12, dpi = 300)