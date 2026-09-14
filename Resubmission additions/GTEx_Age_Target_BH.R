# BH adjustment of the existing seven-gene GTEx grouped-age result table.
# Usage: Rscript GTEx_Age_Target_BH.R input.xlsx output.csv [sheet]
# Input: 7 genes x 6 region groups x 3 tests = 126 rows.
# The default is GTEx_Age_SLIT_ROBO_only, followed by historical sheet names.
# An explicit sheet must still contain the complete seven-gene target table.
# This script adjusts stored raw P values; it does not refit any statistical test.
# BH7 is within region/test, BH35 pools the five regions within test, and BH42
# also includes the pooled All block within test. Test types remain separate.
# These are alternative families; the reporting family needs scientific justification.
# The full ECM-glycoprotein age screen is not represented by this target table.
# Dependencies: readxl and base R. Original row order and source fields are retained.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) stop("Usage: Rscript GTEx_Age_Target_BH.R input.xlsx output.csv [sheet]")
if (file.exists(args[2])) stop("Output already exists; choose a new output path.")
suppressPackageStartupMessages(library(readxl))
sheets <- excel_sheets(args[1])
candidates <- c("GTEx_Age_SLIT_ROBO_only", "GTEx_Age_Raw", "GTEx_Age_Matrisome_Adjusted")
sheet <- if (length(args) == 3L) args[3] else candidates[candidates %in% sheets][1]
if (is.na(sheet) || !sheet %in% sheets) stop("No matching target sheet; supply its name as the third argument.")
data <- as.data.frame(read_excel(args[1], sheet = sheet, .name_repair = "minimal"),
                      check.names = FALSE)
if (nrow(data) != 126L) stop("Expected the 126-row SLIT/ROBO target table, not the full matrisome age table.")
# Retain historical source-state text when supplied; do not invent it when absent.
if ("Source adjustment state" %in% names(data)) {
  if ("Adjustment state" %in% names(data)) stop("Use only one source adjustment-state header.")
  names(data)[names(data) == "Source adjustment state"] <- "Adjustment state"
}
fields <- c("Region", "Gene", "Contrast", "Test", "Raw p", "Older minus younger estimate",
            "Se difference", "Ci lower", "Ci upper", "Effect size", "Direction")
if ("Adjustment state" %in% names(data)) fields <- c(fields, "Adjustment state")
if (anyDuplicated(names(data)) || !all(fields %in% names(data))) stop("Missing or duplicated source headers.")
# Select by header, so an existing derived BH column cannot shift the raw inputs.
data <- data[fields]
regions <- c("All", "CER", "FRO", "HCP", "HTL", "SN")
genes <- c("SLIT1", "SLIT2", "SLIT3", "ROBO1", "ROBO2", "ROBO3", "ROBO4")
tests <- c("Student", "Welch", "Mann-Whitney")
keys <- data[c("Region", "Gene", "Test")]
if (nrow(data) != 126L || anyNA(keys) || anyDuplicated(keys) ||
    !setequal(data$Region, regions) || !setequal(data$Gene, genes) || !setequal(data$Test, tests)) {
  stop("Expected all 126 unique gene/region/test combinations.")
}
p <- data[["Raw p"]]
if (!is.numeric(p) || any(!is.finite(p)) || any(p < 0 | p > 1)) stop("Raw P values must be finite and between zero and one.")
q7 <- q35 <- q42 <- rep(NA_real_, nrow(data))
for (test in tests) {
  same_test <- data$Test == test
  regional <- same_test & data$Region != "All"
  stopifnot(sum(same_test) == 42L, sum(regional) == 35L)
  q42[same_test] <- p.adjust(p[same_test], method = "BH", n = 42L)
  q35[regional] <- p.adjust(p[regional], method = "BH", n = 35L)
  for (region in regions) {
    selected <- same_test & data$Region == region
    stopifnot(sum(selected) == 7L)
    q7[selected] <- p.adjust(p[selected], method = "BH", n = 7L)
  }
}
names(data)[names(data) == "Adjustment state"] <- "Source adjustment state"
data[["BH7 within region and test"]] <- q7
data[["BH35 five regions within test"]] <- q35
data[["BH42 including All within test"]] <- q42
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
write.csv(data, args[2], row.names = FALSE, na = "")
cat("126 source rows; 357 adjusted values; 21 All rows have no BH35 value.\n")
