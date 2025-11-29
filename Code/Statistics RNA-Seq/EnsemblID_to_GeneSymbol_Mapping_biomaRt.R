# Install and load necessary packages if not already installed
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

# Load the necessary libraries
library(biomaRt)
library(readxl)
library(data.table)
library(openxlsx)

# Set the file path
file_path <- "All_genes_normalized_age&tisues_unpaired_5_sites_v1.xlsx"

# Define a function to convert Excel to CSV if needed
convert_to_csv <- function(file_path) {
  csv_file <- sub("\\.xlsx$", ".csv", file_path)
  if (grepl("\\.xlsx$", file_path)) {
    message("Converting Excel file to CSV...")
    data <- read_excel(file_path)
    fwrite(data, csv_file)
    message("Conversion complete: ", csv_file)
  }
  return(csv_file)
}

# Convert to CSV if necessary
csv_file_path <- convert_to_csv(file_path)

# Read the CSV file using fread for faster performance
data <- fread(csv_file_path)

# Extract ENSEMBL gene IDs from the first column
ensembl_gene_ids <- data[[1]]

# Remove version numbers from ENSEMBL IDs (e.g., ENSG00000141510.17 -> ENSG00000141510)
ensembl_gene_ids_clean <- gsub("\\..*$", "", ensembl_gene_ids)

# Connect to the ENSEMBL database using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve gene names corresponding to the ENSEMBL gene IDs
genes <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  filters = 'ensembl_gene_id',
  values = ensembl_gene_ids_clean,
  mart = ensembl
)

# Create a mapping of ENSEMBL IDs to gene names
gene_name_map <- setNames(genes$external_gene_name, genes$ensembl_gene_id)

# Add the gene names as a new column to the original data frame
data$Gene_Name <- gene_name_map[ensembl_gene_ids_clean]

# Save the updated data frame to an Excel file
output_file <- "All_genes_with_gene_names.xlsx"
write.xlsx(data, output_file, row.names = FALSE)
message("Data saved to: ", output_file)