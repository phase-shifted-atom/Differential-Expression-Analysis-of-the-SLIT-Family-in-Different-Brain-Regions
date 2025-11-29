# Load necessary libraries
library(DESeq2)
library(readr)
library(dplyr)

# Function to replace hyphens with dots
replace_hyphens <- function(x) {
  gsub(pattern = "-", replacement = ".", x = x, fixed = TRUE)
}

# File paths
sample_attributes_file <- "GTEx_Analysis_v8_Annotations_SampleAttributes&PhenotypesDS.txt"
gene_reads_file <- "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"

# Check if files exist
if (!file.exists(sample_attributes_file)) {
  stop("Sample attributes file not found.")
}

if (!file.exists(gene_reads_file)) {
  stop("Gene reads file not found.")
}

# Load the sample attributes and phenotypes data
sample_data <- read.delim(sample_attributes_file, stringsAsFactors = FALSE)

# Select relevant columns
sample_attributes <- select(
  sample_data,
  SUBJID, SEX, AGE, DTHHRDY, SMATSSCR, SMCENTER, SAMPID, SMTS, SMTSD, SMAFRZE,
  SMNABTCH, SMGEBTCH, SMRIN, SMTSISCH, FIVEPAIRED
)

# Filter to keep only brain samples that are RNA-Seq
sample_attributes_braindata <- sample_attributes %>%
  filter(SMTS == "Brain" & SMAFRZE == "RNASEQ")

# Include the additional brain regions
brain_regions <- c(
  "Brain - Cerebellar Hemisphere",
  "Brain - Frontal Cortex (BA9)",
  "Brain - Hippocampus",
  "Brain - Hypothalamus",
  "Brain - Substantia nigra"
)

sample_attributes_braindata_sites <- sample_attributes_braindata %>%
  filter(SMTSD %in% brain_regions)

# Replace hyphens with dots in relevant columns and convert to factors
sample_attributes_braindata_sites$SAMPID <- replace_hyphens(sample_attributes_braindata_sites$SAMPID)
sample_attributes_braindata_sites$SMTSD <- factor(replace_hyphens(sample_attributes_braindata_sites$SMTSD))
sample_attributes_braindata_sites$SMNABTCH <- factor(replace_hyphens(sample_attributes_braindata_sites$SMNABTCH))
sample_attributes_braindata_sites$SMGEBTCH <- factor(replace_hyphens(sample_attributes_braindata_sites$SMGEBTCH))
sample_attributes_braindata_sites$SMCENTER <- factor(replace_hyphens(sample_attributes_braindata_sites$SMCENTER))

# Scale SMRIN and SMTSISCH, and add them as new columns
sample_attributes_braindata_sites$SMRIN_scaled <- scale(as.numeric(sample_attributes_braindata_sites$SMRIN))
sample_attributes_braindata_sites$SMTSISCH_scaled <- scale(as.numeric(sample_attributes_braindata_sites$SMTSISCH))

# Convert SUBJID to a factor
sample_attributes_braindata_sites$SUBJID <- factor(sample_attributes_braindata_sites$SUBJID)

# Categorize AGE into age groups
sample_attributes_braindata_sites$AGE_group <- ifelse(sample_attributes_braindata_sites$AGE %in% c("20-29", "30-39", "40-49", "50-59"), "younger", "older")

# Convert AGE_group to a factor
sample_attributes_braindata_sites$AGE_group <- factor(sample_attributes_braindata_sites$AGE_group)

# Perform PCA on batch-related variables (SMNABTCH, SMGEBTCH, SMCENTER)
# Convert categorical variables to numeric factors for PCA
sample_attributes_braindata_sites$SMNABTCH_num <- as.numeric(factor(sample_attributes_braindata_sites$SMNABTCH))
sample_attributes_braindata_sites$SMGEBTCH_num <- as.numeric(factor(sample_attributes_braindata_sites$SMGEBTCH))
sample_attributes_braindata_sites$SMCENTER_num <- as.numeric(factor(sample_attributes_braindata_sites$SMCENTER))

# Perform PCA on the batch-related variables
batch_data <- data.frame(SMNABTCH_num = sample_attributes_braindata_sites$SMNABTCH_num,
                         SMGEBTCH_num = sample_attributes_braindata_sites$SMGEBTCH_num,
                         SMCENTER_num = sample_attributes_braindata_sites$SMCENTER_num)

pca_batch <- prcomp(batch_data, scale. = TRUE)

# View the summary of explained variance
pca_summary <- summary(pca_batch)
message("Variance explained by PC1: ", pca_summary$importance[2, 1])
message("Variance explained by PC2: ", pca_summary$importance[2, 2])
message("Variance explained by PC3: ", pca_summary$importance[2, 3])

# Optionally select components based on explained variance
cumulative_variance <- cumsum(pca_batch$sdev^2 / sum(pca_batch$sdev^2))
num_components <- which(cumulative_variance >= 0.95)[1]
message(paste("Number of principal components selected:", num_components))

# Extract the first few principal components (PC1, PC2, PC3)
pca_data <- as.data.frame(pca_batch$x)
sample_attributes_braindata_sites$PC1_batch <- pca_data$PC1
sample_attributes_braindata_sites$PC2_batch <- pca_data$PC2
sample_attributes_braindata_sites$PC3_batch <- pca_data$PC3

# Load the gene reads data
GTEx_Analysis_gene_reads <- read_table(
  gene_reads_file,
  col_names = TRUE,
  skip = 2,
  col_types = cols(
    Name = col_character(),
    Description = col_character(),
    .default = col_double()  # Ensure numeric data is loaded
  )
)

# Replace hyphens with dots in gene reads data column names
colnames(GTEx_Analysis_gene_reads) <- replace_hyphens(colnames(GTEx_Analysis_gene_reads))

# Create an index for matching samples
sample_columns <- colnames(GTEx_Analysis_gene_reads)[-c(1, 2)]
index <- sample_columns %in% sample_attributes_braindata_sites$SAMPID

# Subset the gene reads data to keep only matching samples
GTEx_Analysis_gene_reads_sites <- as.data.frame(GTEx_Analysis_gene_reads[, c(1, 2, which(index) + 2)])
rownames(GTEx_Analysis_gene_reads_sites) <- GTEx_Analysis_gene_reads_sites$Name

# Extract count data and convert to numeric matrix
count_data <- as.matrix(GTEx_Analysis_gene_reads_sites[, -c(1, 2)])
rownames(count_data) <- GTEx_Analysis_gene_reads_sites$Name

# Verify that count_data is numeric
if (mode(count_data) != "numeric") {
  stop("count_data is not numeric.")
}

# Ensure that sample_info rows match the columns of count_data
samples <- colnames(count_data)
sample_info <- sample_attributes_braindata_sites %>%
  filter(SAMPID %in% samples) %>%
  arrange(match(SAMPID, samples))
rownames(sample_info) <- sample_info$SAMPID

# Verify that samples match
if (!all(samples == rownames(sample_info))) {
  stop("Mismatch between count data columns and sample info rows.")
}

# Convert variables in design formula to factors and clean up levels
sample_info$SEX <- factor(sample_info$SEX)
sample_info$DTHHRDY <- factor(sample_info$DTHHRDY)
sample_info$SMTSD <- factor(sample_info$SMTSD)
levels(sample_info$SMTSD) <- gsub("[^A-Za-z0-9_.]", "_", levels(sample_info$SMTSD))

# Create DESeqDataSet including the principal components from PCA instead of the original batch variables
ds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = ~ SEX + DTHHRDY + PC1_batch + PC2_batch + PC3_batch + SMRIN_scaled + SMTSISCH_scaled + AGE_group * SMTSD
)

# Filter out genes with no reads
keep_genes <- rowSums(counts(ds)) > 0
ds <- ds[keep_genes, ]

message(paste("Number of genes kept after filtering:", sum(keep_genes)))

# Normalize the data
ds <- estimateSizeFactors(ds)
normalized_counts <- counts(ds, normalized = TRUE)

# Write the normalized counts to a file
output_file <- paste0("Normalized_age&tisues_unpaired_5_sites_v1", ".txt")
write.table(
  normalized_counts,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

message("Normalized counts written to file:", output_file)