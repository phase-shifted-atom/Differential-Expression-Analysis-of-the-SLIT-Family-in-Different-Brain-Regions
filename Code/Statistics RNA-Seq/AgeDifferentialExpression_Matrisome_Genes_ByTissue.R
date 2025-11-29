# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stats)
library(effsize)
library(rstatix)  # This package is used to compute the rank-biserial correlation
library(openxlsx)  # This package allows saving in XLSX format

# Load the Excel file and get the sheet names
file_path <- "R_Matrisome_normalized_age&tisues_unpaired_5_sites_v1.xlsx"
sheet_names <- excel_sheets(file_path)

# Initialize an empty list to store results for all sheets
all_results_list <- list()

# Loop through each sheet in the Excel file
for (sheet in sheet_names) {
  # Load the data from the current sheet
  data <- read_excel(file_path, sheet = sheet)
  
  # Separate age groups based on the renamed column "Age Bracket"
  group2 <- data %>% filter(`Age Bracket` == "20-59")
  group1 <- data %>% filter(`Age Bracket` == "60-79")
  
  # Initialize an empty data frame to store results for the current sheet
  sheet_results <- data.frame(Gene = character(),
                              Test = character(),
                              P_Value = numeric(),
                              Location_Parameter = numeric(),
                              SE_Difference = numeric(),
                              CI_Lower = numeric(),
                              CI_Upper = numeric(),
                              Effect_Size = numeric(),
                              Log2_FC = numeric(),
                              stringsAsFactors = FALSE)
  
  # Loop through each gene (column) and perform tests
  for (gene in names(data)[-1]) {  # Exclude the Age Bracket column
    # Extract gene expression values for both age groups
    expr_group1 <- group1[[gene]]
    expr_group2 <- group2[[gene]]
    
    # Perform the Student's t-test
    student_test <- t.test(expr_group1, expr_group2, var.equal = TRUE)
    # Perform the Welch's t-test
    welch_test <- t.test(expr_group1, expr_group2, var.equal = FALSE)
    # Perform the Mann-Whitney U test with Hodges-Lehmann estimator and confidence interval
    mann_whitney_test <- tryCatch(
      wilcox.test(expr_group1, expr_group2, conf.int = TRUE, exact = FALSE),
      warning = function(w) {
        list(p.value = NA, estimate = NA, conf.int = c(NA, NA))  # Return NA for tied data
      },
      error = function(e) {
        list(p.value = NA, estimate = NA, conf.int = c(NA, NA))  # Return NA for errors
      }
    )
    
    # Compute means for log2 fold change
    mean_group1 <- mean(expr_group1, na.rm = TRUE)
    mean_group2 <- mean(expr_group2, na.rm = TRUE)
    
    # Compute Glass's Delta for Welch's t-test (use group2 as control group)
    mean_diff <- mean_group1 - mean_group2
    sd_control <- sd(expr_group2, na.rm = TRUE)
    glass_delta_effect_size <- ifelse(sd_control == 0, NA, mean_diff / sd_control)  # Handle division by zero
    
    # Compute rank-biserial correlation for Mann-Whitney test using rstatix
    combined_data <- data.frame(group = c(rep("20-59", length(expr_group1)), rep("60-79", length(expr_group2))),
                                value = c(expr_group1, expr_group2))
    rank_biserial_effect_size <- wilcox_effsize(combined_data, value ~ group)$effsize
    
    # Compute log2 fold change
    # If one of the means is zero or NA, log2 fold change will be NA
    if (!is.na(mean_group1) && !is.na(mean_group2) && mean_group2 != 0) {
      log2_FC <- log2(mean_group1 / mean_group2)
    } else {
      log2_FC <- NA
    }
    
    # Store Student's t-test results
    sheet_results <- rbind(sheet_results, data.frame(Gene = gene,
                                                     Test = "Student",
                                                     P_Value = student_test$p.value,
                                                     Location_Parameter = student_test$estimate[1] - student_test$estimate[2],
                                                     SE_Difference = student_test$stderr,
                                                     CI_Lower = student_test$conf.int[1],
                                                     CI_Upper = student_test$conf.int[2],
                                                     Effect_Size = cohen.d(expr_group1, expr_group2)$estimate,
                                                     Log2_FC = log2_FC))
    
    # Store Welch's t-test results with Glass's Delta
    sheet_results <- rbind(sheet_results, data.frame(Gene = gene,
                                                     Test = "Welch",
                                                     P_Value = welch_test$p.value,
                                                     Location_Parameter = welch_test$estimate[1] - welch_test$estimate[2],
                                                     SE_Difference = welch_test$stderr,
                                                     CI_Lower = welch_test$conf.int[1],
                                                     CI_Upper = welch_test$conf.int[2],
                                                     Effect_Size = glass_delta_effect_size,
                                                     Log2_FC = log2_FC))
    
    # Store Mann-Whitney test results with Hodges-Lehmann estimator and rank-biserial correlation
    sheet_results <- rbind(sheet_results, data.frame(Gene = gene,
                                                     Test = "Mann-Whitney",
                                                     P_Value = mann_whitney_test$p.value,
                                                     Location_Parameter = mann_whitney_test$estimate,  # Hodges-Lehmann estimator
                                                     SE_Difference = NA,  # Mann-Whitney does not compute SE
                                                     CI_Lower = mann_whitney_test$conf.int[1],  # 95% CI lower bound
                                                     CI_Upper = mann_whitney_test$conf.int[2],  # 95% CI upper bound
                                                     Effect_Size = rank_biserial_effect_size,
                                                     Log2_FC = log2_FC))
  }
  
  # Adjust p-values within each test type for the current sheet
  sheet_results <- sheet_results %>%
    group_by(Test) %>%
    mutate(Adjusted_P_Value = p.adjust(P_Value, method = "BH")) %>%
    ungroup()
  
  # Add the sheet name to the results
  sheet_results$Sheet <- sheet
  
  # Reorder columns
  sheet_results <- sheet_results[, c("Sheet", "Gene", "Test", "P_Value", "Adjusted_P_Value",
                                     "Location_Parameter", "SE_Difference", "CI_Lower", "CI_Upper", "Effect_Size", "Log2_FC")]
  
  # Append the results to the list
  all_results_list[[sheet]] <- sheet_results
}

# Combine all sheets' results into one data frame
all_results <- do.call(rbind, all_results_list)

# Save the combined results to an XLSX file in the project folder
write.xlsx(all_results, "GTEX_Age_Analysis_Results_All_Sheets_with_BH_Log2FC.xlsx", row.names = FALSE)