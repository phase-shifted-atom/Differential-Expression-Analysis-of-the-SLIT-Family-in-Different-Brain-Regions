# ============================================================
# R Script: Parallelized Analysis of Transposed Gene Expression Data
# With Added log2FC Column
# ============================================================

# -------------------------------
# 1. Load Necessary Libraries
# -------------------------------
required_packages <- c("readxl", "dplyr", "tidyr", "stats", 
                       "effsize", "rstatix", "openxlsx",
                       "parallel", "foreach", "doParallel")  # PARALLEL

installed_packages <- rownames(installed.packages())
for(pkg in required_packages){
  if(!pkg %in% installed_packages){
    install.packages(pkg, dependencies = TRUE)
  }
}

library(readxl)
library(dplyr)
library(tidyr)
library(stats)
library(effsize)
library(rstatix)
library(openxlsx)
library(parallel)   # PARALLEL
library(foreach)    # PARALLEL
library(doParallel) # PARALLEL

# -------------------------------
# 2. Define File Paths and Parameters
# -------------------------------
data_file_path <- "R_Matrisome_normalized_age&tisues_unpaired_5_sites_v1"  # Input Excel file
output_file    <- "GTEX_Age_Analysis_Results_All_Sheets_with_BH.xlsx"  # Output results file

# Define your two group labels (case-insensitive)
group_labels <- c("CTL", "LGG")  # e.g., c("GroupA", "GroupB")

# -------------------------------
# 3. Retrieve Sheet Names
# -------------------------------
sheet_names <- excel_sheets(data_file_path)
cat("Sheets found in the Excel file:", paste(sheet_names, collapse = ", "), "\n\n")

# -------------------------------
# 4. Initialize Lists to Store Results and Diagnostics
# -------------------------------
all_results_list <- list()

diagnostic_info <- data.frame(
  Sheet = character(),
  Total_Genes = integer(),
  Columns_group1 = integer(),  # Renamed for general usage
  Columns_group2 = integer(),  # Renamed for general usage
  Unique_Group_Labels = character(),
  stringsAsFactors = FALSE
)

# -------------------------------
# 5. Initialize Parallel Backend
# -------------------------------
## PARALLEL: Detect and register cores
n_cores <- detectCores() - 1  # or a fixed number
cat("Using", n_cores, "cores for parallel execution.\n")
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# -------------------------------
# 6. Loop Through Each Sheet
# -------------------------------
for (sheet in sheet_names) {
  cat("Processing sheet:", sheet, "\n")
  
  # -------------------------------
  # 6.a. Read the Sheet Data
  # -------------------------------
  data_initial <- tryCatch(
    read_excel(path = data_file_path, sheet = sheet, col_names = FALSE),
    error = function(e) {
      warning(paste("Sheet:", sheet, "- Error reading the sheet. Skipping this sheet."))
      return(NULL)
    }
  )
  
  if (is.null(data_initial)) {
    next  # Skip to the next sheet if there's an error
  }
  
  # -------------------------------
  # 6.b. Extract Group Labels from First Row
  # -------------------------------
  first_row <- data_initial[1, -1] 
  group_labels_in_sheet <- trimws(gsub("[^ -~]", "", as.character(first_row[1, ])))
  
  unique_labels <- unique(group_labels_in_sheet)
  cat("  Unique Group Labels in sheet:", paste(unique_labels, collapse = ", "), "\n")
  
  # -------------------------------
  # 6.c. Identify Columns for Each Group
  # -------------------------------
  group_columns <- list()
  
  for (lbl in group_labels) {
    cols <- which(tolower(group_labels_in_sheet) == tolower(lbl))
    if (length(cols) == 0) {
      warning(paste("Sheet:", sheet, "- No columns found for group label:", lbl))
    } else {
      group_columns[[lbl]] <- cols
      cat("  Identified", length(cols), "columns for group label:", lbl, "\n")
    }
  }
  
  # Validate that both groups were found
  missing_groups <- setdiff(group_labels, names(group_columns))
  if (length(missing_groups) > 0) {
    warning(
      paste(
        "Sheet:", sheet, 
        "- The following group labels were not found:", 
        paste(missing_groups, collapse = ", "), 
        "- Skipping this sheet."
      )
    )
    # Record partial diagnostics
    diagnostic_info <- rbind(
      diagnostic_info,
      data.frame(
        Sheet = sheet,
        Total_Genes = NA,
        Columns_group1 = ifelse(!is.null(group_columns[[group_labels[1]]]),
                                length(group_columns[[group_labels[1]]]), 0),
        Columns_group2 = ifelse(!is.null(group_columns[[group_labels[2]]]),
                                length(group_columns[[group_labels[2]]]), 0),
        Unique_Group_Labels = paste(unique_labels, collapse = ", "),
        stringsAsFactors = FALSE
      )
    )
    next
  }
  
  # Assign group1 and group2
  group1_label <- group_labels[1]
  group2_label <- group_labels[2]
  
  group1_cols <- group_columns[[group1_label]]
  group2_cols <- group_columns[[group2_label]]
  
  # -------------------------------
  # 6.d. Extract Gene Names and Expression Data
  # -------------------------------
  gene_names <- as.character(data_initial[-1, 1][[1]])
  total_genes <- length(gene_names)
  cat("  Total number of genes:", total_genes, "\n")
  
  expression_data <- data_initial[-1, -1]
  expression_matrix <- as.matrix(expression_data)
  
  # -------------------------------
  # 6.e. Parallel Loop Through Each Gene
  # -------------------------------
  # We'll create a single data frame (sheet_results) from combining
  # the small data frames returned by each iteration.
  
  # PARALLEL: We'll do the big for-loop over genes using foreach.
  # Each iteration returns 1-3 rows for Student, Welch, Mann-Whitney.
  
  sheet_results <- foreach(
    i = seq_len(total_genes), 
    .combine = rbind,          # row-bind results
    .packages = c("effsize","rstatix","stats") 
    # 'dplyr' not strictly necessary in loop if we do data.frame manipulation
  ) %dopar% {
    
    gene <- gene_names[i]
    
    # Extract expression for both groups
    expr_group1 <- as.numeric(expression_matrix[i, group1_cols])
    expr_group2 <- as.numeric(expression_matrix[i, group2_cols])
    
    valid_group1 <- sum(!is.na(expr_group1))
    valid_group2 <- sum(!is.na(expr_group2))
    if (valid_group1 < 2 | valid_group2 < 2) {
      # Return empty if not enough data
      return(NULL)
    }
    
    # Initialize placeholders
    student_p  <- NA
    welch_p    <- NA
    mann_whitney_p <- NA
    
    student_estimate      <- NA
    welch_estimate        <- NA
    mann_whitney_estimate <- NA
    
    student_se <- NA
    student_ci_lower <- NA
    student_ci_upper <- NA
    
    welch_se <- NA
    welch_ci_lower <- NA
    welch_ci_upper <- NA
    
    cohen_d <- NA
    glass_delta <- NA
    rank_biserial <- NA
    
    # ---------- Student's t-test ----------
    student_test <- tryCatch(
      t.test(expr_group1, expr_group2, var.equal = TRUE),
      error = function(e) NULL
    )
    if (!is.null(student_test)) {
      student_p         <- student_test$p.value
      # t.test normally returns estimate c(mean of x, mean of y)
      student_estimate  <- student_test$estimate[1] - student_test$estimate[2]
      student_se        <- student_test$stderr
      student_ci_lower  <- student_test$conf.int[1]
      student_ci_upper  <- student_test$conf.int[2]
      
      cohen_d_result <- tryCatch(
        cohen.d(expr_group1, expr_group2)$estimate,
        error = function(e) NA
      )
      cohen_d <- cohen_d_result
    }
    
    # ---------- Welch's t-test ----------
    welch_test <- tryCatch(
      t.test(expr_group1, expr_group2, var.equal = FALSE),
      error = function(e) NULL
    )
    
    mean_group1 <- NA
    mean_group2 <- NA
    
    if (!is.null(welch_test)) {
      welch_p         <- welch_test$p.value
      welch_estimate  <- welch_test$estimate[1] - welch_test$estimate[2]
      welch_se        <- welch_test$stderr
      welch_ci_lower  <- welch_test$conf.int[1]
      welch_ci_upper  <- welch_test$conf.int[2]
      
      mean_group1 <- mean(expr_group1, na.rm = TRUE)
      mean_group2 <- mean(expr_group2, na.rm = TRUE)
      sd_control  <- sd(expr_group2, na.rm = TRUE)
      
      if (sd_control != 0) {
        glass_delta <- (mean_group1 - mean_group2) / sd_control
      } else {
        glass_delta <- NA
      }
    } else {
      # If Welch test errored, compute these anyway for log2_FC
      mean_group1 <- mean(expr_group1, na.rm = TRUE)
      mean_group2 <- mean(expr_group2, na.rm = TRUE)
    }
    
    # ---------- Mann-Whitney U Test ----------
    mann_whitney_test <- tryCatch(
      wilcox.test(expr_group1, expr_group2, conf.int = TRUE, exact = FALSE),
      warning = function(w) {
        list(p.value = NA, estimate = NA, conf.int = c(NA, NA))
      },
      error = function(e) {
        list(p.value = NA, estimate = NA, conf.int = c(NA, NA))
      }
    )
    mann_whitney_p         <- mann_whitney_test$p.value
    mann_whitney_estimate  <- mann_whitney_test$estimate
    mann_whitney_ci_lower  <- mann_whitney_test$conf.int[1]
    mann_whitney_ci_upper  <- mann_whitney_test$conf.int[2]
    
    # Rank-Biserial Correlation
    # We'll treat group1_label as "group1" and group2_label as "group2"
    # for ordering. (Accessing them is possible as global variables in parallel:
    # doParallel tries to export objects automatically, but you might need
    # .export = c("group1_label","group2_label") in foreach if there's an error.)
    combined_data <- data.frame(
      group = factor(c(
        rep(group1_label, length(expr_group1)),
        rep(group2_label, length(expr_group2))
      )),
      value = c(expr_group1, expr_group2)
    )
    
    rank_biserial_result <- tryCatch(
      wilcox_effsize(combined_data, value ~ group)$effsize,
      error = function(e) NA
    )
    rank_biserial <- rank_biserial_result
    
    # ---------- Compute log2_FC ----------
    if (!is.na(mean_group1) && !is.na(mean_group2) && mean_group2 != 0) {
      log2_FC <- log2(mean_group1 / mean_group2)
    } else {
      log2_FC <- NA
    }
    
    # We might have up to 3 rows to return: Student, Welch, Mann-Whitney
    # We'll gather them in a small list of data frames and rbind at the end.
    result_rows <- list()
    
    # Student's t-test row
    if (!is.na(student_p)) {
      result_rows[[length(result_rows)+1]] <- data.frame(
        Sheet              = sheet,
        Gene               = gene,
        Test               = "Student",
        P_Value            = student_p,
        Adjusted_P_Value   = NA,  # Fill later
        Location_Parameter = student_estimate,
        SE_Difference      = student_se,
        CI_Lower           = student_ci_lower,
        CI_Upper           = student_ci_upper,
        Effect_Size        = cohen_d,
        Log2_FC            = log2_FC,
        stringsAsFactors   = FALSE
      )
    }
    
    # Welch's t-test row
    if (!is.na(welch_p)) {
      result_rows[[length(result_rows)+1]] <- data.frame(
        Sheet              = sheet,
        Gene               = gene,
        Test               = "Welch",
        P_Value            = welch_p,
        Adjusted_P_Value   = NA,
        Location_Parameter = welch_estimate,
        SE_Difference      = welch_se,
        CI_Lower           = welch_ci_lower,
        CI_Upper           = welch_ci_upper,
        Effect_Size        = glass_delta,
        Log2_FC            = log2_FC,
        stringsAsFactors   = FALSE
      )
    }
    
    # Mann-Whitney row (always appended, because we always do that test)
    result_rows[[length(result_rows)+1]] <- data.frame(
      Sheet              = sheet,
      Gene               = gene,
      Test               = "Mann-Whitney",
      P_Value            = mann_whitney_p,
      Adjusted_P_Value   = NA,
      Location_Parameter = mann_whitney_estimate,
      SE_Difference      = NA,
      CI_Lower           = mann_whitney_ci_lower,
      CI_Upper           = mann_whitney_ci_upper,
      Effect_Size        = rank_biserial,
      Log2_FC            = log2_FC,
      stringsAsFactors   = FALSE
    )
    
    # Combine the small list into one data frame for this gene
    if (length(result_rows) > 0) {
      do.call(rbind, result_rows)
    } else {
      NULL
    }
  }
  
  cat("  Completed processing genes for sheet:", sheet, 
      " --> Result rows:", nrow(sheet_results), "\n")
  
  # -------------------------------
  # 6.f. Adjust P-Values Using BH
  # -------------------------------
  if (nrow(sheet_results) > 0) {
    sheet_results <- sheet_results %>%
      group_by(Test) %>%
      mutate(Adjusted_P_Value = p.adjust(P_Value, method = "BH")) %>%
      ungroup()
  }
  
  # -------------------------------
  # 6.g. Append Sheet Results
  # -------------------------------
  all_results_list[[sheet]] <- sheet_results
  
  # -------------------------------
  # 6.h. Record Diagnostic Info
  # -------------------------------
  diagnostic_info <- rbind(
    diagnostic_info,
    data.frame(
      Sheet = sheet,
      Total_Genes = total_genes,
      Columns_group1 = length(group1_cols),
      Columns_group2 = length(group2_cols),
      Unique_Group_Labels = paste(unique_labels, collapse = ", "),
      stringsAsFactors = FALSE
    )
  )
  
  cat("  Results appended for sheet:", sheet, "\n\n")
}  # end loop over sheets

# -------------------------------
# 7. Combine All Sheets' Results
# -------------------------------
cat("Combining results from all sheets...\n")
all_results <- bind_rows(all_results_list)
cat("Combination complete. Total results:", nrow(all_results), "rows.\n\n")

# -------------------------------
# 8. Save the Combined Results
# -------------------------------
cat("Saving results to Excel file:", output_file, "\n")
write.xlsx(all_results, output_file, rowNames = FALSE)
cat("Results successfully saved to", output_file, "\n\n")

# -------------------------------
# 9. Save Diagnostic Information
# -------------------------------
diagnostic_file <- "Age_Group_Diagnostic_Report.xlsx"
write.xlsx(diagnostic_info, diagnostic_file, rowNames = FALSE)
cat("Diagnostic information saved to", diagnostic_file, "\n")

# -------------------------------
# 10. Stop the Parallel Cluster
# -------------------------------
## PARALLEL: Clean up
stopCluster(cl)
cat("Parallel cluster stopped.\n")

# ============================
# End of Script
# ============================