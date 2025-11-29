# ---------------------------------------------------------------------------
# Script: GTEx_Age_Lollipop_Combined_MatrisomeAndSites.R
# Purpose:
#   Generate two coordinated lollipop plots summarising age-related
#   differential expression of hippocampal matrisome genes:
#     (1) Genes grouped by Matrisome type (left panel).
#     (2) Genes grouped by brain site (right panel), with SLIT genes
#         highlighted across sites.
#
# Input:
#   - "GTEx Wilcoxon Glycoproteins.xlsx"
#       Table of significant hippocampal matrisome genes with at least:
#         * gene               : gene symbol
#         * Matrisome_type     : Matrisome class/category
#         * Fishers_combined_p : pooled p-value (age effect, Wilcoxon-based)
#   - "Bar Plot Significant Glycoproteins.xlsx"
#       Table of significant glycoprotein genes across brain sites, with:
#         * Gene_Symbol       : gene symbol
#         * Brain_Site        : brain region code (HCP, FRO, HTL, CER)
#         * Adjusted_P_Value  : multiple-comparison adjusted p-value
#
# Output:
#   - "Combined_LollipopCharts.png"
#       Multi-panel figure combining:
#         * left panel: -log10(Fisher's combined p) by gene and Matrisome type
#         * right panel: -log10(adjusted p) by gene and brain site
#
# Notes:
#   - Both panels:
#       * Use consistent axis label sizes (axis_label_size) for readability.
#       * Highlight SLIT1/2/3 labels in red and add a red outline around the
#         corresponding points, as these genes are a focus of the manuscript.
#   - This script assumes upstream differential expression / Wilcoxon tests
#     have already been run and that the input tables contain only the genes
#     of interest for visualisation.
# ---------------------------------------------------------------------------

# Combined script: build and align the two lollipop plots described above

# Load necessary libraries
# - readxl      : import Excel spreadsheets with the differential results
# - ggplot2     : build the lollipop plots
# - dplyr       : data wrangling (mutate, arrange, group_by)
# - scales      : transformations if needed (e.g., -log10)
# - ggtext      : render HTML-formatted axis labels (colored / italic text)
# - patchwork   : combine the two ggplot objects into a single figure
library(readxl)
library(ggplot2)
library(dplyr)
library(scales)
library(ggtext)     # for element_markdown()
library(patchwork)  # for combining plots

# ---- Common settings ----
# Use a single axis_label_size so that both panels share identical axis
# text sizing, which makes the combined figure visually coherent.
axis_label_size <- 8    # use the same size for both x and y axis labels

# ---- Plot 1: Genes by Matrisome Type (left panel) ----
# This panel shows, for hippocampal matrisome genes, the strength of the
# age effect (–log10(Fisher's combined p)) coloured by Matrisome type.

# Read the matrisome summary table (hippocampus-only Wilcoxon results)
data1 <- read_excel("GTEx Wilcoxon Glycoproteins.xlsx")

# Standardise column names: replace spaces with underscores to simplify
# column access in dplyr/ggplot pipelines.
names(data1) <- gsub(" ", "_", names(data1))

# Create:
#   - gene_cluster : unique identifier combining gene symbol and Matrisome type
#   - negLogP      : –log10(Fisher's combined p-value) for plotting on x-axis
data1 <- data1 %>%
  mutate(gene_cluster = paste(gene, Matrisome_type, sep = "_"),
         negLogP      = -log10(Fishers_combined_p))

# Order Matrisome types by their maximum significance across genes
# (types with more strongly associated genes appear first in the legend).
Matrisome_type_order <- data1 %>%
  group_by(Matrisome_type) %>%
  summarise(max_negLogP = max(negLogP, na.rm = TRUE)) %>%
  arrange(desc(max_negLogP)) %>%
  pull(Matrisome_type)

# Within each Matrisome type:
#   1. Arrange genes by p-value (most significant first).
#   2. Set gene_cluster as a factor and reverse the levels so that the
#      most significant genes appear at the top of the y-axis.
data1 <- data1 %>%
  arrange(match(Matrisome_type, Matrisome_type_order),
          Fishers_combined_p) %>%
  mutate(gene_cluster = factor(gene_cluster,
                               levels = rev(unique(gene_cluster))))

# Determine a common left x-axis limit for this panel
# (used as the baseline from which segments or points extend).
x_min1 <- floor(min(data1$negLogP, na.rm = TRUE))

# Define consistent colour palettes for Matrisome types:
# - matrisome_fills   : softer colours for point interiors
# - matrisome_strokes : more saturated colours for point borders and legend
matrisome_fills <- c(
  "ECM-affiliated Proteins" = "#A1AFB6",  # azul-ardósia claro
  "ECM Glycoproteins"       = "#C187D0",  # lilás claro
  "ECM Regulators"          = "#CBA48B",  # tan claro
  "Proteoglycans"           = "#EA82A8",  # framboesa clara
  "Secreted Factors"        = "#AF9D96",  # taupe claro
  "Collagens"               = "#AEA5C5"   # lavanda acinzentada
)
matrisome_strokes <- c(
  "ECM-affiliated Proteins" = "#546E7A",  # azul-cinza médio
  "ECM Glycoproteins"       = "#8E24AA",  # roxo vívido
  "ECM Regulators"          = "#A05A2C",  # cobre
  "Proteoglycans"           = "#D81B60",  # carmim
  "Secreted Factors"        = "#6D4C41",  # mocha
  "Collagens"               = "#6B5B95"   # ametista
)

# Prepare y-axis labels: strip the Matrisome_type suffix so that only
# the gene symbols appear, then mark SLIT genes in red via HTML.
base_labels1    <- sub("_.*", "", levels(data1$gene_cluster))
colored_labels1 <- ifelse(
  base_labels1 %in% c("SLIT1","SLIT2","SLIT3"),
  paste0("<span style='color:red;'>", base_labels1, "</span>"),
  base_labels1
)

# Build the left lollipop plot:
#   - each point corresponds to a gene × Matrisome type combination
#   - SLIT1/2/3 get an extra red outline for emphasis
p_left <- ggplot(data1, aes(x = negLogP, y = gene_cluster, color = Matrisome_type)) +
  # Main points: filled circles coloured by Matrisome class
  geom_point(aes(fill = Matrisome_type, color = Matrisome_type), shape = 21, size = 3, stroke = 0.9) +
  # Extra outline for SLIT1/2/3 genes (focus of the manuscript)
  geom_point(data = filter(data1, gene %in% c("SLIT1","SLIT2","SLIT3")),
             aes(x = negLogP, y = gene_cluster),
             shape = 21, fill = NA, color = "red", size = 5, stroke = 1) +
  # Do not extend the x-axis further left than the minimum –log10(p)
  scale_fill_manual(values = matrisome_fills) +
  scale_color_manual(values = matrisome_strokes) +
  scale_x_continuous(limits = c(x_min1, NA)) +
  # Apply HTML-formatted y-axis labels (ggtext::element_markdown will render them)
  scale_y_discrete(labels = colored_labels1) +
  labs(
    x    = "-log10(Fisher's Combined p-value)",
    y    = "Gene",
    fill = "Matrisome Type"
  ) +
  # Use a minimal theme and apply the shared axis_label_size to both axes
  # so that this panel visually matches the right-hand panel.
  theme_minimal() +
  theme(
    axis.text.y   = element_markdown(size = axis_label_size),
    axis.text.x   = element_text(size = axis_label_size),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(ncol = 3), color = "none")

# ---- Plot 2: Genes by Brain Site with SLIT Highlight (right panel) ----
# This panel shows, for significant glycoprotein genes, the strength of the
# association at each brain site (–log10 adjusted p-value), with SLIT genes
# highlighted across sites.

# Specify the display order of brain sites on the legend and y-axis labels
# (ensuring consistency with other figures and the main text).
brain_site_order <- c("HCP", "FRO", "HTL", "CER")  # display order

# Read the cross-site glycoproteins summary table
data2 <- read_excel("Bar Plot Significant Glycoproteins.xlsx")

# Standardise column names as in Plot 1 (spaces → underscores)
names(data2) <- gsub(" ", "_", names(data2))

# Construct:
#   - Gene_Brain_Site : unique identifier for each gene × site combination
#   - Highlight       : logical flag for SLIT1/2/3 (for red outline)
#   - Brain_Site      : ordered factor to enforce the desired site order
# Then:
#   - arrange rows by Brain_Site and Adjusted_P_Value (most significant first)
#   - turn Gene_Brain_Site into a factor and reverse levels to put the most
#     significant entries toward the top of the y-axis.
data2 <- data2 %>%
  mutate(
    Gene_Brain_Site = paste(Gene_Symbol, Brain_Site, sep = "_"),
    Highlight       = Gene_Symbol %in% c("SLIT1","SLIT2","SLIT3"),
    Brain_Site      = factor(Brain_Site, levels = brain_site_order)
  ) %>%
  arrange(Brain_Site, Adjusted_P_Value) %>%
  mutate(
    Gene_Brain_Site = factor(Gene_Brain_Site, levels = rev(Gene_Brain_Site))
  )

# Colour palette for brain sites (used for points and legend keys)
brain_site_colors <- c(
  "HCP" = "#ff9300",
  "HTL" = "#ff8ad8",
  "FRO" = "#7a81ff",
  "CER" = "#83cceb"
)

# Build the right lollipop-style plot:
#   - horizontal segments from 0 to –log10(adjusted p)
#   - coloured points by brain site
#   - vertical dotted line at the conventional 0.05 threshold
p_right <- ggplot(data2, aes(y = Gene_Brain_Site, x = -log10(Adjusted_P_Value), color = Brain_Site)) +
  # Baseline segments: encode the magnitude of –log10(adjusted p)
  geom_segment(aes(yend = Gene_Brain_Site, x = 0, xend = -log10(Adjusted_P_Value)),
               color = "gray70") +
  # Main points: coloured by brain site
  geom_point(size = 3) +
  # Significance reference line at p = 0.05 (–log10(0.05))
  geom_vline(xintercept = -log10(0.05), linetype = "dotted") +
  # Custom y-axis labels: extract gene symbols and render SLIT1/2/3 in red
  scale_color_manual(values = brain_site_colors) +
  scale_y_discrete(labels = function(x) {
    genes <- sub("_.*", "", x)
    ifelse(genes %in% c("SLIT1","SLIT2","SLIT3"),
           paste0("<span style='color:red;'>", genes, "</span>"),
           genes)
  }) +
  # Start the x-axis at 1 (–log10(0.1)) so that extremely small bars on the
  # left are not visually compressed; add a small expansion on the right.
  scale_x_continuous(limits = c(1, NA), expand = c(0, 0.1)) +
  labs(
    x     = "-log10 Adjusted p-value",
    y     = NULL,
    color = "Brain Site"
  ) +
  # Use the same minimal theme and axis_label_size as in the left panel
  # so that both panels look visually integrated when combined.
  theme_minimal() +
  theme(
    axis.text.y   = element_markdown(size = axis_label_size),
    axis.text.x   = element_text(size = axis_label_size),
    legend.position = "right"
  ) +
  # Add a red outline around SLIT1/2/3 points in all brain sites
  geom_point(data = filter(data2, Highlight),
             aes(x = -log10(Adjusted_P_Value), y = Gene_Brain_Site),
             shape = 21, fill = NA, color = "red", size = 4, stroke = 1)


# ---- Combine & save ----

# Combine the two panels side-by-side
# (right panel first, then left, to match the manuscript figure layout).
combined <- p_right + p_left + plot_layout(ncol = 2)

# Display the combined figure in the current plotting device
print(combined)

# Save a high-resolution PNG of the combined figure
# Width/height tuned to maintain legible axis labels and legends.
ggsave("Combined_LollipopCharts.png",
       plot  = combined,
       width = 10,
       height = 14,
       dpi    = 300)