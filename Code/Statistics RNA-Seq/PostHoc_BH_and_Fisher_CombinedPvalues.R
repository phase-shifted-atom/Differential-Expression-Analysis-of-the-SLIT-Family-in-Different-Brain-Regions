# Load necessary libraries
library(readxl)
library(dplyr)
library(openxlsx)
library(foreach)
library(doParallel)

# Function to install and load missing packages
install_and_load <- function(packages) {
  installed_packages <- rownames(installed.packages())
  for (pkg in packages) {
    if (!(pkg %in% installed_packages)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

# Install and load required packages
required_packages <- c("readxl", "dplyr", "openxlsx", "foreach", "doParallel")
install_and_load(required_packages)

# Define file path
file_path <- "R_Matrisome_normalized_only_tissues_paired_5_sites_v1_Transposed.xlsx"

# Read data without setting column names
cat("Lendo o arquivo Excel...\n")
data_raw <- read_excel(file_path, col_names = FALSE)
cat("Leitura concluída.\n")

# Extract gene names (first column, excluding the first row)
cat("Extraindo nomes dos genes...\n")
gene_names <- data_raw %>%
  slice(-1) %>%         # Remove the first row (header of brain sites)
  pull(1) %>%           # Extract the first column as a vector
  as.character()        # Convert to character
cat("Número de genes extraídos:", length(gene_names), "\n")

# Extract brain site labels (first row, excluding the first column)
cat("Extraindo rótulos dos sítios cerebrais...\n")
brain_sites <- data_raw %>%
  slice(1) %>%           # Select the first row
  select(-1) %>%         # Remove the first column
  unlist() %>%           # Convert to vector
  as.character()         # Ensure they are characters
cat("Número de sítios cerebrais extraídos:", length(brain_sites), "\n")

# Extract gene expression data (removing first row and first column)
cat("Extraindo dados de expressão gênica...\n")
expression_data <- data_raw %>%
  slice(-1) %>%          # Remove the first row
  select(-1) %>%         # Remove the first column
  as.matrix()            # Convert to matrix
mode(expression_data) <- "numeric"
cat("Dimensões dos dados de expressão (linhas x colunas):",
    dim(expression_data)[1], "x", dim(expression_data)[2], "\n")

# Transpose the expression matrix
cat("Transpondo a matriz de expressão...\n")
expression_transposed <- as.data.frame(t(expression_data))
cat("Dimensões após transposição (linhas x colunas):",
    dim(expression_transposed)[1], "x", dim(expression_transposed)[2], "\n")

# Assign gene names as column names
if (length(gene_names) != ncol(expression_transposed)) {
  stop("O número de nomes de genes não corresponde ao número de colunas no dataframe transposto.")
}
colnames(expression_transposed) <- gene_names

# Check for duplicated gene names and make them unique
if (any(duplicated(colnames(expression_transposed)))) {
  cat("Nomes de genes duplicados encontrados. Tornando-os únicos.\n")
  colnames(expression_transposed) <- make.unique(colnames(expression_transposed))
} else {
  cat("Todos os nomes de genes são únicos.\n")
}

# Add 'Brain Site' column with the extracted labels
cat("Adicionando a coluna 'Brain Site'...\n")
expression_transposed <- expression_transposed %>%
  mutate(`Brain Site` = brain_sites) %>%
  select(`Brain Site`, everything())

# Ensure 'Brain Site' is a factor
expression_transposed <- expression_transposed %>%
  mutate(`Brain Site` = as.factor(`Brain Site`))

# Check the final structure of transposed data
cat("Estrutura final de expression_transposed:\n")
print(str(expression_transposed))
cat("Primeiras linhas de expression_transposed:\n")
print(head(expression_transposed))

# Define all reference brain sites and their respective comparison sites
reference_sites <- list(
  "Hippocampus"             = c("Frontal Cortex (BA9)", "Cerebellar Hemisphere", "Substantia nigra", "Hypothalamus"),
  "Frontal Cortex (BA9)"    = c("Hippocampus", "Cerebellar Hemisphere", "Substantia nigra", "Hypothalamus"),
  "Hypothalamus"            = c("Hippocampus", "Frontal Cortex (BA9)", "Cerebellar Hemisphere", "Substantia nigra"),
  "Substantia nigra"        = c("Hippocampus", "Frontal Cortex (BA9)", "Cerebellar Hemisphere", "Hypothalamus"),
  "Cerebellar Hemisphere"   = c("Hippocampus", "Frontal Cortex (BA9)", "Substantia nigra", "Hypothalamus")
)

# Set up parallel processing
num_cores <- parallel::detectCores() - 1  # Reserve one core for the system
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("Processamento paralelo configurado com", num_cores, "núcleos.\n")

# Start Wilcoxon (Mann–Whitney) tests in parallel
cat("Iniciando os testes de Wilcoxon em paralelo...\n")

all_results_df <- foreach(reference = names(reference_sites), .combine = 'rbind', .packages = c('dplyr')) %dopar% {
  
  comparison_sites <- reference_sites[[reference]]
  temp_results <- vector("list", length = length(comparison_sites) * length(gene_names))
  index <- 1
  gene_count <- 0
  total_genes <- length(gene_names)
  
  for (gene in gene_names) {
    for (comparison in comparison_sites) {
      
      # Subset data for current brain site and reference
      site_data <- expression_transposed %>%
        filter(`Brain Site` %in% c(comparison, reference)) %>%
        select(`Brain Site`, !!sym(gene)) %>%
        mutate(`Brain Site` = factor(`Brain Site`, levels = c(comparison, reference)))
      
      # Check if there is enough data for the test
      if (length(unique(site_data$`Brain Site`)) < 2) {
        p_val <- NA
      } else {
        # Run one-sided Wilcoxon test (less alternative)
        test_result <- tryCatch(
          wilcox.test(
            site_data[[gene]] ~ site_data$`Brain Site`,
            alternative = "less",
            exact = FALSE,
            correct = FALSE
          ),
          error = function(e) NA
        )
        p_val <- ifelse(is.list(test_result), test_result$p.value, NA)
      }
      
      temp_results[[index]] <- data.frame(
        Reference_Site    = reference,
        Comparison        = paste(comparison, "vs", reference),
        Gene              = gene,
        p_value           = p_val,
        stringsAsFactors  = FALSE
      )
      index <- index + 1
    }
    
    gene_count <- gene_count + 1
    if (gene_count %% 100 == 0) {
      cat(paste0("Processado ", gene_count, " genes para referência ", reference, " de ", total_genes, " genes...\n"))
    }
  }
  
  cat(paste0("Concluído o processamento para referência: ", reference, "\n"))
  do.call(rbind, temp_results)
}

# Stop parallel cluster
stopCluster(cl)
cat("Processamento paralelo concluído.\n")

# Check results structure
cat("Número total de resultados:", nrow(all_results_df), "\n")

# Apply Benjamini-Hochberg correction per reference site
cat("Aplicando correção de Benjamini-Hochberg...\n")
all_results_df <- all_results_df %>%
  group_by(Reference_Site) %>%
  mutate(p_value_adjusted = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Compute Fisher combined p-value for each Gene and Reference_Site
cat("Calculando valores-p combinados de Fisher...\n")
fisher_pvals <- all_results_df %>%
  group_by(Reference_Site, Gene) %>%
  summarise(Fishers_combined_p = {
    valid_pvals <- p_value_adjusted[!is.na(p_value_adjusted)]
    if (length(valid_pvals) > 0) {
      chi_sq <- -2 * sum(log(valid_pvals))
      pchisq(chi_sq, df = 2 * length(valid_pvals), lower.tail = FALSE)
    } else {
      NA
    }
  }, .groups = 'drop')

# Merge Fisher p-values back
all_results_df <- all_results_df %>%
  left_join(fisher_pvals, by = c("Reference_Site", "Gene"))

# Reorder columns
all_results_df <- all_results_df %>%
  select(Reference_Site, Comparison, Gene, p_value, p_value_adjusted, Fishers_combined_p)

# Create output directory
output_dir <- "Wilcoxon_Test_Results_Global_BH_Per_Reference"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Save results per reference site
cat("Salvando os resultados em arquivos Excel separados...\n")
for (ref in unique(all_results_df$Reference_Site)) {
  ref_results <- all_results_df %>%
    filter(Reference_Site == ref) %>%
    arrange(Gene, Comparison)
  
  file_name   <- paste0(gsub("[() ]", "_", ref), "_Wilcoxon_test_results_with_BH_and_Fisher.xlsx")
  output_path <- file.path(output_dir, file_name)
  
  tryCatch({
    write.xlsx(ref_results, output_path, rowNames = FALSE)
    cat("Resultados salvos para", ref, "em", output_path, "\n")
  }, error = function(e) {
    cat("Erro ao salvar para", ref, ":", e$message, "\n")
  })
}

cat("Todos os resultados foram salvos com sucesso.\n")