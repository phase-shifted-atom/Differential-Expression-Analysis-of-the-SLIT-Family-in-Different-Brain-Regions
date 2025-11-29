# Repository: Statistical analysis code for RNA-seq, immunohistochemistry (IHC) and graph production

## Title of the study

**Differential Expression Analysis of the SLIT Family in Different Brain Regions**

### Authors, affiliations, and contact

- Authors blinded for review.

### Citation and links

- Journal: Proceedings of the Royal Society B: Biological Sciences

---

## Short summary of the study

The study combines bulk RNA-sequencing of post-mortem human brain tissue from a local biorepository (n = 63 samples across eight regions) with external validation in the GTEx v8 cohort (n = 962 samples across five regions) to chart the expression landscape of the 239-member core matrisome in the ageing human brain. Hippocampus emerges as a major ECM outlier, with 114 matrisome genes—dominated by ECM glycoproteins—up-regulated relative to seven other regions. Within this ECM landscape, SLIT1–3 and their receptor ROBO2 rank among the most strongly enriched signatures and show robust co-expression, identifying a hippocampus-centric, age-labile SLIT–ROBO2 module.

Age-stratified analyses in GTEx (younger < 60 years vs older ≥ 60 years) reveal a concerted late-life decline of SLIT–ROBO2 expression, accompanied by downregulation of the PAK–LIMK–cofilin arm of the pathway that stabilises F-actin. This suggests that attenuation of SLIT–ROBO2 signalling weakens synaptic scaffolding and contributes to hippocampal vulnerability during cognitive ageing. Immunohistochemistry (IHC) and exploratory western blotting further support region-specific gradients and subcellular compartmentalisation of SLIT1–3 and ROBO2 proteins across hippocampus, cortex, cerebellum, hypothalamus and substantia nigra.

This repository provides:

- **RNA-seq analysis code** for preprocessing, normalisation, differential expression, pathway-level summaries and figure generation for both the local biorepository and GTEx cohorts.
- **IHC analysis code** (ImageJ/Fiji macros) for image preprocessing, colour deconvolution, ROI quantification, statistical analyses and visualisation.
- **Graph production scripts** R files for generating plots, heatmaps and related figures from the processed RNA-seq outputs.
- **Header files and example outputs** that document the structure of all input tables and reproduce the numerical source data underlying Figures 1–4 and Supplementary Figures 1–3 of the paper.

---

## Data availability, requests, and header files

### Data availability and how to request data

The RNA-seq data analysed in this study can be found in GEO under submission Series GSE311762.

- **Local biorepository RNA-seq and protein data**
  - Local transcriptomic data comprise 63 post-mortem brain samples across eight regions (CER, HCP, HTL, SN, TEMP, OCC, FRO, PAR), sequenced as 3′ mRNA libraries and processed to gene-level counts. Demographic and clinical metadata include sex, age, post-mortem interval, cause of death, comorbidities and lifestyle factors. IHC and western blot data were generated on matched regions from independent donors.
  - These data are governed by institutional ethics approvals and cannot be shared directly via this repository. Researchers interested in accessing the underlying raw or processed data should contact the corresponding author, who can provide guidance on the necessary ethics approvals and data-use agreements with the institutional biobank.

- **GTEx v8 RNA-seq data**
  - External validation analyses use GTEx v8 bulk RNA-seq data from five brain regions (CER, FRO, HCP, HTL and SN), including both a paired subset for regional comparisons and a larger age-stratified cohort (< 60 vs ≥ 60 years).
  - Access to GTEx RNA-seq data should be requested via the GTEx/dbGaP portal, following their controlled-access procedures.

- **IHC images and western blots**
  - Brightfield IHC images (SLIT1–3, ROBO2) and western blot data were generated from paraffin-embedded or frozen tissue from three donors per region, after quality-control exclusions. High-resolution images, raw microscope files and full western blot scans are governed by biobank and ethics regulations and are not deposited here.

Researchers aiming to reproduce the analyses should:

1. Obtain access to GTEx v8 count-level data for the relevant brain regions via the GTEx/dbGaP procedures.
2. Access the local biorepository RNA-seq in GEO.
3. Use the **header files and example tables** in this repository to confirm that their locally obtained datasets match the expected input structure before running the scripts.

This repository is designed so that, once equivalent raw or count-level data and IHC images are obtained from the original sources, all analyses and figures/tables in the paper can be regenerated using the scripts provided here.

---

### Header files: structure and units

To document the expected input format, we provide header files in `data_headers/`. These files contain either the header row alone or the header row plus a few anonymised/example rows. Below we describe the main headers and corresponding units.

#### RNA-seq header files

Header/example files in `data_headers/` document the expected input structure for each R script. All header files are provided as comma-separated values (`.csv`) and include either only the header row or the header plus a small number of anonymised example rows.

- **`AgeDifferentialExpression_Matrisome_Genes_ByTissue Headers.csv`**  
  Structure of the input matrices used by `AgeDifferentialExpression_Matrisome_Genes_ByTissue.R` for age-stratified matrisome analyses.  
  - `Age Bracket`: Age group label for each sample (e.g. `20–59`, `60–79`).  
    *Unit:* years, grouped into categorical brackets.  
  - One column per gene (e.g. `SEMA3F`, `WNT16`, `PLXND1`, …): normalised expression values for matrisome genes in the tissue under analysis, with column names given as HGNC gene symbols.  
    *Unit:* DESeq2 size-factor–normalised counts (dimensionless, relative expression scale).  
  - Each row corresponds to a single sample within a given tissue.

- **`TwoGroupDE_AllGenes_Parallelized Headers.csv`**  
  Structure of the input matrices used by `TwoGroupDE_AllGenes_Parallelized.R` for general two-group comparisons. The layout matches the file above and can be reused whenever a binary grouping variable is available.  
  - `Age Bracket` (or another grouping variable, depending on the analysis): group label for each sample.  
    *Unit:* categorical.  
  - One column per gene (HGNC gene symbol): normalised expression values for all genes included in the analysis.  
    *Unit:* same expression scale used for the corresponding dataset (typically DESeq2-normalised counts; dimensionless).  
  - Each row corresponds to one sample.

- **`EnsemblID_to_GeneSymbol_Mapping_biomaRt Headers.csv`**  
  Structure of gene-level expression tables accepted by `EnsemblID_to_GeneSymbol_Mapping_biomaRt.R` when adding HGNC symbols to Ensembl IDs.  
  - `ensembl_gene_id`: Ensembl gene identifier without version suffix (e.g. `ENSG00000223972`).  
    *Unit:* identifier (no unit).  
  - One column per sample (e.g. `GTEX.1192X.0011.R10a.SM.DO941`, `GTEX.11DXW.0011.R1a.SM.DNZZD`, …): expression values or other gene-level statistics that will be carried through unchanged when gene symbols are appended.  
    *Unit:* depends on the upstream pipeline (for the provided example: DESeq2 normalised counts, dimensionless).  
  - Each row corresponds to a single gene.

- **`Limma_DifferentialExpression_BrainSites_SVA_SOM data Headers.csv`**  
  Structure of the gene-by-sample count matrix used by `Limma_DifferentialExpression_BrainSites_SVA_SOM.R`.  
  - `Unnamed: 0`: Gene identifier (gene symbol as used in the DGEList object, e.g. `A1BG`, `A2M`).  
    *Unit:* identifier (no unit).  
  - `width`: Gene length (as returned by the annotation step used to build the DGEList; e.g. exon span in base pairs).  
    *Unit:* base pairs.  
  - One column per sample (e.g. `10_S10`, `11_S11`, …, `63_S63`): raw gene-level counts per sample.  
    *Unit:* integer read counts per gene.  
  - Each row corresponds to a single gene.

- **`Limma_DifferentialExpression_BrainSites_SVA_SOM pData Headers.csv`**  
  Structure of the sample-level metadata (phenotype data frame) used together with the matrix above.  
  - `Unnamed: 0`: Sample identifier matching the column names in the count matrix (e.g. `X10_S10`).  
    *Unit:* identifier (no unit).  
  - `Site`: Brain region label for each sample (e.g. `HCP`, `FRO`, `HTL`, `OCC`, `TEMP`, `CER`, `SN`, `PAR`).  
    *Unit:* categorical.  
  - `Case`: Donor/case identifier.  
    *Unit:* identifier (no unit).  
  - `Gender`: Biological sex of the donor (`M`, `F`).  
    *Unit:* categorical.  
  - `Age`: Age of the donor at death.  
    *Unit:* years.  
  - `DM2`: Indicator of type 2 diabetes status (`Y` / `N`).  
    *Unit:* categorical (yes/no).

- **`PostHoc_BH_and_Fisher_CombinedPvalues Headers.csv`**  
  Structure of wide tables used by `PostHoc_BH_and_Fisher_CombinedPvalues.R` when aggregating gene-level measures across brain regions or related contrasts.  
  - `Brain Site`: Gene identifier (for the example file, an HGNC gene symbol such as `SEMA3F`, `WNT16`, `PLXND1`, …).  
    *Unit:* identifier (no unit).  
  - Subsequent columns correspond to region-specific values for each donor, grouped by brain region and, where applicable, suffixed by a numeric index to distinguish repeated measures (e.g. `Frontal Cortex (BA9)`, `Cerebellar Hemisphere`, `Hippocampus`, `Substantia nigra`, `Hypothalamus`, followed by `Frontal Cortex (BA9).1`, `Cerebellar Hemisphere.1`, …).  
    *Unit:* gene-level summary statistics (e.g. normalised expression values or other continuous measures) on the scale produced by the preceding analysis step; dimensionless.  
  - Each row corresponds to one gene.

#### IHC header files

All IHC macros run on .TIFF images, not on spreadsheets.

Quantitative measurements exported from ImageJ/Fiji are subsequently analysed in R; header templates for these tables are provided in `data_headers/`:

- **`IHC_ROI_measurements_header.tsv`**  
  ROI-level measurements exported from ImageJ/Fiji after running the colour deconvolution and measurement macros. Typical columns include:  
  - `Image` or `File`: image filename.  
    *Unit:* identifier (no unit).  
  - `ROI` or `Label`: region-of-interest identifier.  
    *Unit:* identifier (no unit).  
  - `Area`: area of the ROI in square micrometres.  
    *Unit:* µm².  
  - `Mean` or `Mean_intensity`: mean pixel intensity within the ROI for the stain of interest.  
    *Unit:* arbitrary intensity units (a.u.).  
  - Additional columns may include standard ImageJ measurement outputs (e.g. `Min`, `Max`, `StdDev`) depending on which options were selected in the *Set Measurements* dialog.  
    *Unit:* as defined by ImageJ (dimensionless or a.u.).

- **`IHC_SummaryPerCase_header.tsv`**  
  Per-case and per-region summaries derived from the ROI-level tables and used as input for IHC statistical analyses in R. Typical columns include:  
  - `CaseID`: donor/case identifier.  
    *Unit:* identifier (no unit).  
  - `Region`: brain region or anatomical subregion (e.g. hippocampus, cortex, cerebellum).  
    *Unit:* categorical.  
  - `Marker`: stain or antibody (e.g. `SLIT1`, `SLIT2`, `SLIT3`, `ROBO2`).  
    *Unit:* categorical.  
  - `Metric`: summary measure name (e.g. mean positive area, H-score, mean intensity).  
    *Unit:* categorical label.  
  - `Value`: numerical summary statistic for the given case/region/marker/metric.  
    *Unit:* depends on the metric (e.g. µm² for area-based metrics, a.u. for intensity-based metrics, dimensionless for H-scores).  
  - Optional columns can store additional covariates used in the analyses (e.g. age group, diagnosis) as categorical variables.

#### Graph header files

Header/example files in `data_headers/` also define the structure of the tables used as direct numerical source data for the lollipop plots and heatmaps generated by the scripts in `code/`. Each file is a `.csv` containing the minimum set of columns needed to recreate the corresponding figure panel (plus, where appropriate, a few anonymised example rows).

- **`GTEx_Age_Lollipop_ByMatrisomeType Headers.csv`**  
  Structure of the per-gene statistics tables used by `GTEx_Age_Lollipop_ByMatrisomeType.R` for age-stratified matrisome lollipop plots.  
  - `gene`: HGNC gene symbol (e.g. `COL25A1`, `RSPO3`, `PCSK5`).  
    *Unit:* identifier (no unit).  
  - `Matrisome_group`: High-level matrisome grouping (e.g. `Core matrisome`, `Matrisome-associated`).  
    *Unit:* categorical.  
  - `Matrisome_type`: Matrisome category (e.g. `Collagens`, `ECM Glycoproteins`, `ECM Regulators`, `ECM-affiliated Proteins`, `Proteoglycans`).  
    *Unit:* categorical.  
  - `p_value`: Nominal p-value for the primary age comparison in the tissue of interest.  
    *Unit:* probability (0–1, dimensionless).  
  - `p_value_adjusted`: Multiple-testing–corrected p-value (e.g. Benjamini–Hochberg FDR) for the same comparison.  
    *Unit:* probability (0–1, dimensionless).  
  - `Fishers_combined_p`: Fisher’s combined p-value summarising evidence across related tests for each gene.  
    *Unit:* probability (0–1, dimensionless; smaller values indicate stronger aggregate evidence).  
  - Each row corresponds to a single gene.

- **`GTEx_Age_Lollipop_Combined_MatrisomeAndSites A Headers.csv`**  
  Structure of the summary tables used by `GTEx_Age_Lollipop_Combined_MatrisomeAndSites.R` to generate the panel that combines matrisome type with an across-site summary.  
  - `gene`: HGNC gene symbol.  
    *Unit:* identifier (no unit).  
  - `Matrisome_type`: Matrisome category for each gene (e.g. `ECM Glycoproteins`, `Proteoglycans`).  
    *Unit:* categorical.  
  - `Fishers_combined_p`: Fisher’s combined p-value across site-specific or related age tests, used to rank genes within each matrisome type.  
    *Unit:* probability (0–1, dimensionless).  
  - Each row corresponds to one gene and provides the statistics used to position that gene in the combined matrisome-type lollipop plot.

- **`GTEx_Age_Lollipop_Combined_MatrisomeAndSites B Headers.csv`**  
  Structure of the per-site summary tables used by `GTEx_Age_Lollipop_Combined_MatrisomeAndSites.R` to generate the panel that displays site-wise adjusted p-values.  
  - `Brain Site`: Brain region label (e.g. `HCP`, `HTL`, `FRO`, `CER`, `SN`).  
    *Unit:* categorical.  
  - `Gene Symbol`: HGNC gene symbol for the matrisome gene of interest (e.g. `NELL1`, `SLIT1`, `RELN`).  
    *Unit:* identifier (no unit).  
  - `Adjusted_P_Value`: Multiple-testing–corrected p-value (e.g. Benjamini–Hochberg FDR) for the age comparison within the indicated brain site and gene.  
    *Unit:* probability (0–1, dimensionless).  
  - Each row corresponds to one gene–site pair.

- **`Matrisome_Zscore_Heatmap_BrainSites Headers.csv`**  
  Structure of the expression matrices used by `Matrisome_Zscore_Heatmap_BrainSites.R` to generate heatmaps of matrisome expression across brain sites.  
  - `Group`: Matrisome group label for each gene (e.g. `ECM Glycoproteins`, `ECM Regulators`, `Proteoglycans`, `Collagens`, `Secreted Factors`).  
    *Unit:* categorical.  
  - `Gene`: HGNC gene symbol.  
    *Unit:* identifier (no unit).  
  - Columns beginning with brain-region prefixes such as `CER`, `SUB`, `HTL`, `HCP`, `TEMP`, `OCC`, `FRO`, `PAR`, `SN` (with optional numeric suffixes `.1`, `.2`, …) contain expression values per donor and region, for example `CER`, `CER.1`, …, `CER.7` for cerebellar samples and `PAR`, `PAR.1`, …, `PAR.15` for parietal cortex samples.  
    *Unit:* expression values on the same scale as the upstream pipeline (e.g. normalised counts; dimensionless). These columns are subsequently centred and scaled (Z-scored) by the plotting script when building the heatmap.  
  - Each row corresponds to a single gene; each column after `Gene` corresponds to a region-specific replicate used to build the heatmap.

- **`ROBO_Pathway_Age_ZscoreHeatmap Headers.csv`**  
  Structure of the age-stratified expression tables used by `ROBO_Pathway_Age_ZscoreHeatmap.R` for SLIT–ROBO pathway heatmaps.  
  - `Group`: Pathway branch label (e.g. `Ena/VASP`, `LIMK`) or functional submodule within the SLIT–ROBO signalling cascade.  
    *Unit:* categorical.  
  - `Age Bracket`: Column storing the HGNC gene symbol within each pathway branch (e.g. `EVL`, `VASP`, `CFL1`, `DSTN`).  
    *Unit:* identifier (no unit).  
  - Columns whose names start with `20-59` (e.g. `20-59`, `20-59.1`, …) contain per-sample expression values for individuals in the younger age group (20–59 years).  
    *Unit:* expression values on the same scale as the upstream pipeline (e.g. DESeq2-normalised counts; dimensionless).  
  - Columns whose names start with `60-79` (e.g. `60-79`, `60-79.1`, …) contain per-sample expression values for individuals in the older age group (60–79 years).  
    *Unit:* expression values on the same scale as the upstream pipeline (e.g. DESeq2-normalised counts; dimensionless).  
  - The heatmap script applies robust scaling/Z-scoring to these per-sample columns before plotting.  
  - Each row corresponds to one pathway gene within a specific pathway branch.

---

## Software environment: code and package versions

All analyses were performed using **R** and **ImageJ/Fiji**, with additional graphical work in **GraphPad Prism**.

- **R:** `R version 4.4.1` with Bioconductor (e.g. `Bioconductor 3.20`).  
- **ImageJ/Fiji:** Fiji distribution (ImageJ 1.54g) with the built-in *Colour Deconvolution* plugin.  
- **GraphPad Prism:** GraphPad Prism 9 (used for Welch estimation plots).

- **Key R packages used across scripts (not exhaustive):**
  - Core statistical / RNA-seq packages:  
    `DESeq2`, `edgeR`, `limma`, `sva`, `Biobase`, `genefilter`
  - Annotation and import:  
    `org.Hs.eg.db`, `AnnotationDbi`, `tximport`, `biomaRt`, `readr`, `readxl`, `data.table`, `openxlsx`
  - Data handling and I/O:  
    `dplyr`, `tidyr`
  - Effect size and statistics helpers (RNA-seq and IHC):  
    `effsize`, `rstatix`, `gtools`
  - Visualisation:  
    `ggplot2`, `pheatmap`, `randomcoloR`, `ggbiplot`
  - Dimension reduction and modelling (RNA-seq):  
    `supraHex`, `GSVA`, `fastICA` or equivalent ICA implementation
  - Parallelisation:  
    `parallel`, `foreach`, `doParallel`

---

## Repository structure and file overview

├── code/
│   ├── DESeq2_Normalization_Brain_Paired5Sites.R
│   ├── DESeq2_Normalization_Brain_AgeByTissue_Unpaired.R
│   ├── AgeDifferentialExpression_Matrisome_Genes_ByTissue.R
│   ├── TwoGroupDE_AllGenes_Parallelized.R
│   ├── PostHoc_BH_and_Fisher_CombinedPvalues.R
│   ├── EnsemblID_to_GeneSymbol_Mapping_biomaRt.R
│   ├── Limma_DifferentialExpression_BrainSites_SVA_SOM.R
│   ├── GTEx_Age_Lollipop_ByMatrisomeType.R
│   ├── GTEx_Age_Lollipop_Combined_MatrisomeAndSites.R
│   ├── Matrisome_Zscore_Heatmap_BrainSites.R
│   ├── ROBO_Pathway_Age_ZscoreHeatmap.R
│   ├── IHC_Background_Removal.ijm
│   ├── ImageJ_BatchCrop_SaveROIs.ijm
│   └── ImageJ_ColorDeconvolution_BatchMeasurements.ijm
├── data_headers/
│   ├── AgeDifferentialExpression_Matrisome_Genes_ByTissue Headers.csv
│   ├── TwoGroupDE_AllGenes_Parallelized Headers.csv
│   ├── EnsemblID_to_GeneSymbol_Mapping_biomaRt Headers.csv
│   ├── Limma_DifferentialExpression_BrainSites_SVA_SOM data Headers.csv
│   ├── Limma_DifferentialExpression_BrainSites_SVA_SOM pData Headers.csv
│   ├── PostHoc_BH_and_Fisher_CombinedPvalues Headers.csv
│   ├── GTEx_Age_Lollipop_ByMatrisomeType Headers.csv
│   ├── GTEx_Age_Lollipop_Combined_MatrisomeAndSites A Headers.csv
│   ├── GTEx_Age_Lollipop_Combined_MatrisomeAndSites B Headers.csv
│   ├── Matrisome_Zscore_Heatmap_BrainSites Headers.csv
│   ├── ROBO_Pathway_Age_ZscoreHeatmap Headers.csv
│   ├── IHC_ROI_measurements_header.tsv
│   └── IHC_SummaryPerCase_header.tsv
└── README.md

### Code files

Below we summarise the purpose of each script. All scripts are **annotated with comments explaining the rationale** for each analysis step. Where necessary, the code has been structured to make the implementation as clear as possible.

#### RNA-seq analysis scripts

- **`code/DESeq2_Normalization_Brain_Paired5Sites.R`**  
  Performs DESeq2 normalisation on gene-level count data for brain samples restricted to subjects with paired samples across the selected five brain sites (GTEx paired design).  
  - **Input:** GTEx (or equivalent) gene read counts and sample attributes.  
  - **Key steps:** filter to brain samples; restrict to subjects with complete sampling in the chosen sites; summarise batch variables via PCA; build a DESeq2 design that includes site and batch PCs; filter genes with no reads; compute normalised counts.  
  - **Output:** Normalised count matrices used as input for downstream paired-site comparisons and for deriving source data for figures/tables involving paired region analyses.

- **`code/DESeq2_Normalization_Brain_AgeByTissue_Unpaired.R`**  
  Performs DESeq2 normalisation on gene-level count data including all selected brain tissues (unpaired), with an age-by-tissue design.  
  - **Input:** GTEx (or equivalent) gene read counts and sample attributes.  
  - **Key steps:** filter to brain samples; define age groups and tissue factors; summarise batch variables (SMNABTCH, SMGEBTCH, SMCENTER, etc.) via PCA; include the resulting PCs in the DESeq2 design; filter low-count genes; compute normalised counts.  
  - **Output:** Normalised count matrices for each brain tissue and age group, used in age-related matrisome/non-matrisome analyses and corresponding figure/table source data.

- **`code/AgeDifferentialExpression_Matrisome_Genes_ByTissue.R`**  
  Performs differential expression analyses between age groups for matrisome genes (and, optionally, non-matrisome genes) across tissues.  
  - **Input:** Excel file(s) containing normalised expression values (one sheet per tissue or dataset).  
  - **Key steps:** for each tissue/sheet, divide samples into age groups; run, per gene, Student’s t-test, Welch’s t-test and Mann–Whitney U test; compute Hodges–Lehmann estimates, effect sizes (e.g. Glass’s Δ, rank-biserial correlation) and log2 fold change; adjust p-values using Benjamini–Hochberg; compile results across tissues.  
  - **Output:** Tables of per-gene statistics by tissue and age group, which serve as source data for age-related figures/tables.

- **`code/TwoGroupDE_AllGenes_Parallelized.R`**  
  General-purpose pipeline for two-group comparisons across all genes (e.g. control vs disease, or any binary grouping).  
  - **Input:** Excel file(s) with expression matrices (one sheet per dataset), with columns labelled to indicate group membership.  
  - **Key steps:** automatically detect columns for each group; for each gene and sheet, perform Student’s t-test, Welch’s t-test and Mann–Whitney U test; compute mean expression per group, log2 fold change and effect sizes; apply Benjamini–Hochberg correction; use parallelisation (`foreach`, `doParallel`) to speed up computation for large gene sets.  
  - **Output:** Combined result tables with per-gene statistics for each dataset, used as source data for figures/tables involving general two-group comparisons.

- **`code/PostHoc_BH_and_Fisher_CombinedPvalues.R`**  
  Post-hoc script for additional multiple-testing control and meta-analysis across comparisons.  
  - **Input:** Results tables where rows correspond to genes and columns to p-values for different site-specific tests or related comparisons.  
  - **Key steps:** reformat tables as needed; apply Benjamini–Hochberg correction within reference sites or groups; compute Fisher’s combined p-values across multiple related tests per gene; organise results by reference; export to Excel.  
  - **Output:** Per-gene adjusted and combined p-values, which can be used to support summary statistics and source data tables for figures/tables that integrate multiple comparisons.

- **`code/EnsemblID_to_GeneSymbol_Mapping_biomaRt.R`**  
  Utility script to add gene symbols to tables based on Ensembl IDs.  
  - **Input:** Table(s) with Ensembl gene IDs (optionally including version suffixes) and associated statistics.  
  - **Key steps:** strip version suffixes; query Ensembl using `biomaRt`; append HGNC gene symbols and other annotation fields; export the annotated table.  
  - **Output:** Annotated result tables with gene symbols, used in final source data for figures and tables.

- **`code/Limma_DifferentialExpression_BrainSites_SVA_SOM.R`**  
  Voom/limma-based differential expression analysis comparing brain sites, with SVA and SOM/heatmap visualisation.  
  - **Input:** Gene-level counts and sample metadata (including brain site labels, sex and potential covariates).  
  - **Key steps:** assemble a DGEList and apply TMM normalisation; filter lowly expressed genes; run voom to obtain log2-CPM and precision weights; use SVA to estimate unmodelled technical variation and build an extended design; fit limma models; construct all pairwise contrasts between brain sites; summarise numbers of differentially expressed genes; generate heatmaps of differentially expressed genes; perform SOM clustering (`supraHex`) to identify expression patterns; run GO/KEGG/MSigDB enrichment; compute GSVA scores and ICA on expression matrices.  
  - **Output:** Differential expression tables, SOM cluster assignments, enrichment results, GSVA scores and ICA-based summaries, which constitute the source data for figures/tables involving site-specific patterns and pathway-level analyses.

#### IHC image analysis and statistical scripts

- **`code/IHC_Background_Removal.ijm`**  
  ImageJ/Fiji macro used to harmonise pale, low-stain tissue images while keeping the background pure white.  
  - **Purpose:** Standardise background and low-stain regions across batches to facilitate consistent thresholding and colour deconvolution, especially in slides with high non-specific staining.  
  - **Key steps:** convert images to HSB space; define tunable ranges for background and low-stain tissue (hue/saturation/brightness); generate masks; optionally shrink masks; adjust saturation/brightness in the low-stain component; optionally run in batch mode over a folder of images.  
  - **Input/Output:** Reads brightfield IHC images (e.g. TIFF) and saves processed images with uniform white background and diminished low-stain visibility.

- **`code/ImageJ_BatchCrop_SaveROIs.ijm`**  
  ImageJ/Fiji macro to batch-crop large TIFF images using a fixed rectangular ROI and save the cropped images.  
  - **Purpose:** Enforce a consistent field-of-view across all IHC images (e.g. focusing on a specific anatomical layer or region) before quantification.  
  - **Key steps:** prompt for input/output folders; allow the user to adjust the fixed ROI coordinates; optionally recurse into subfolders; count eligible files for progress reporting; crop each image to the chosen ROI; save cropped images with a configurable suffix; optionally overwrite or skip existing outputs.  
  - **Input/Output:** Reads full-resolution TIFF images and writes cropped TIFFs suitable for downstream deconvolution and measurement.

- **`code/ImageJ_ColorDeconvolution_BatchMeasurements.ijm`**  
  ImageJ/Fiji macro to perform colour deconvolution and batch measurements on all images in a folder.  
  - **Purpose:** Separate stains (e.g. haematoxylin and DAB) and extract quantitative measurements per image/ROI for use in statistical analyses.  
  - **Key steps:** prompt for an input folder; iterate over eligible images; run colour deconvolution (e.g. built-in *Colour Deconvolution* command) with the relevant stain matrix; generate stain-specific images; apply ROIs (from ROI Manager or pre-saved ROI sets) and call `Measure` on the desired channels; append results to the ImageJ *Results* table for export.  
  - **Input/Output:** Reads cropped IHC images and produces ROI-level measurements (e.g. area, mean intensity) that are exported as tab-delimited tables and used as input to the R IHC statistics scripts.

- **`code/IHC_Statistics_SummarisePerCase_And_GroupComparisons.R`**  
  Main R script for IHC statistical analyses.  
  - **Purpose:** Convert ROI-level measurements into biologically interpretable per-case/per-region metrics and perform the statistical summaries reported in the manuscript.  
  - **Key steps:** import ROI-level tables produced by ImageJ; compute summary statistics per stain/region/case (e.g. mean or median positive area, H-score); merge with clinical/phenotype data (e.g. age group, diagnosis, region); compute descriptive statistics and effect sizes as appropriate; correct for multiple testing (e.g. Benjamini–Hochberg) when inferential analyses are run; export tidy tables ready to be used as figure/table source data.  
  - **Output:** Tidy per-case and per-region summary tables used as numerical source data for IHC figures and tables.

#### Figure generation scripts (RNA-seq)
 
- **`GTEx_Age_Lollipop_ByMatrisomeType.R`**  
  Produces lollipop plots summarising age-related differential expression in hippocampal matrisome genes grouped by matrisome category (e.g. ECM glycoproteins, proteoglycans, collagens). These plots are used for panels that display age-stratified enrichment of ECM glycoproteins in GTEx brain regions.
  - **Input:** Excel table of significant genes with columns such as `gene`, `Matrisome_type`, p-values (e.g. `p_value` or `padj`) and, optionally, log2 fold changes.  
  - **Key steps:** tidy and re-code matrisome categories; create a gene–category identifier; compute –log10(p) statistics; select top genes per category; set factor levels to control ordering; generate a lollipop plot using custom colour palettes and save it to file.  
  - **Output:** Vector and/or raster figure files for the corresponding manuscript figure panels.
 
- **`GTEx_Age_Lollipop_Combined_MatrisomeAndSites.R`**  
  Produces side-by-side lollipop plots combining matrisome-type information and brain-site information in a single figure, as used for multi-region age-related comparisons.  
  - **Input:** As above, plus additional columns indicating brain site or tissue.  
  - **Key steps:** define shared graphical parameters (axis label size, palettes); build two coordinated lollipop plots (matrisome-type and site-wise summaries); harmonise theme elements (fonts, tick labels, legend); use `patchwork` to arrange panels side-by-side.  
  - **Output:** Combined multi-panel figure file and associated source data tables for each panel.
 
- **`Matrisome_Zscore_Heatmap_BrainSites.R`**  
  Generates a heatmap of Z-scored expression values for significant matrisome genes across brain sites, used for visualising core matrisome profiles across local and GTEx regions.
  - **Input:** Excel file with matrisome class, gene name and per-site Z-scores (e.g. `Limma Matrisome HCP Significant Wilcoxon v2.xlsx`).  
  - **Key steps:** extract classes/genes and build an expression matrix; optionally Winsorise extreme values; define ordered factors for site and matrisome group; build annotation data frames and colour palettes; call `pheatmap` with a symmetric colour scale; save the heatmap and, optionally, the underlying ordered matrix.  
  - **Output:** Heatmap image(s) and the ordered, Winsorised Z-score matrix used as figure source data.
 
- **`ROBO_Pathway_Age_ZscoreHeatmap.R`**  
  Re-usable “universal heatmap” script used here for SLIT2–ROBO2 pathway age comparisons (e.g. Figure 4B).
  - **Input:** Excel file of pathway genes with per-group or per-age Z-scores (e.g. `HCP ROBO2 Signaling Significant Genes v6.xlsx`).  
  - **Key steps:** define robust scaling and Winsorisation settings; order rows and columns according to biologically defined groups (e.g. pathway branches, age brackets); perform robust row scaling; optionally detect and down-weight outliers; build annotations and colour maps; generate heatmaps with `pheatmap` and export an Excel workbook with the scaled data.  
  - **Output:** generate heatmaps and an accompanying workbook with the processed matrix, suitable as source data for the relevant figure(s).
 
> **Note:** During development some scripts may have had working names (e.g. `Combined GTEx Expression and Age Fig 1.R`, `Lollipop Plot Grouped by Matrisome Type v4.R`, `Matrisome Heatmap.R`, `ROBO pathway age heatmap.R`). For the final repository, they have been renamed to descriptive names so that the mapping between scripts and figures is explicit.

---

## Instructions for users to run the software

This section outlines a typical workflow to reproduce the analyses after obtaining the relevant data.

### 1. Clone the repository and set up software

- Clone or download this repository.
- Install the required R packages listed above (e.g. using `install.packages` and `BiocManager::install` as appropriate).
- Install Fiji (ImageJ distribution) and verify that the *Colour Deconvolution* plugin is available.

### 2. Prepare data files

- **RNA-seq**  
  - Obtain the relevant count matrices and sample metadata from the original sources (e.g. GTEx v8 via dbGaP, institutional/local biorepository datasets).  
  - Make sure that:  
    - Sample identifiers in the count matrices and metadata match (`SAMPID` / column names).  
    - Tissues, age groups and other covariates are coded consistently with the assumptions in the scripts.  
  - Optionally, adapt the header files in `data_headers/` to confirm that your data conform to the expected structure.

- **IHC**  
  - Obtain the raw or pre-exported IHC images and associated phenotype/clinical data from the corresponding biobank or project.  
  - Ensure that image filenames and case identifiers match the IDs used in the clinical/phenotype tables and in the ROI measurement files.  
  - Use the header files in `data_headers/` as a template to verify that your tables have the expected column names and units.

### 3. Run DESeq2 normalisation (RNA-seq)

- For **paired brain site analyses** (GTEx paired 5-site dataset):  
  - Edit file paths at the top of `code/DESeq2_Normalization_Brain_Paired5Sites.R` to point to your local count and metadata files.  
  - Run the script in R/RStudio. It will:  
    - Filter to brain samples with complete sampling in the chosen sites.  
    - Summarise batch variables via PCA.  
    - Perform DESeq2 normalisation and save normalised count tables.

- For **age-by-tissue analyses (unpaired)**:  
  - Edit file paths in `code/DESeq2_Normalization_Brain_AgeByTissue_Unpaired.R`.  
  - Run the script to obtain normalised counts stratified by tissue and age group.

### 4. Age-related differential expression (RNA-seq)

- Assemble input matrices (e.g. Excel files) with one sheet per tissue, containing the normalised expression values for matrisome/non-matrisome genes.
- Edit file paths and sheet lists at the top of `code/AgeDifferentialExpression_Matrisome_Genes_ByTissue.R` as needed.
- Run the script to obtain per-gene statistics by tissue and age group. The output tables can be saved in `figure_source_data/` and referenced as source data for the corresponding figures/tables.

### 5. General two-group differential expression (RNA-seq)

- For disease vs control or other binary comparisons, organise your data as an Excel file with one sheet per dataset.
- Adjust the group labels and sheet list at the top of `code/TwoGroupDE_AllGenes_Parallelized.R`.
- Run the script to produce per-gene statistics, which can also be saved as figure/table source data.

### 6. Brain site comparisons and clustering (RNA-seq)

- Prepare a count matrix and metadata table for brain site analyses as specified in `data_headers/`.
- Edit file paths at the top of `code/Limma_DifferentialExpression_BrainSites_SVA_SOM.R`.
- Run the script. It will:  
  - Normalise counts using edgeR (TMM) and voom.  
  - Estimate surrogate variables with SVA.  
  - Fit limma models for all pairwise site contrasts.  
  - Generate differential expression tables, heatmaps, SOM-based clusters, enrichment results, GSVA scores and ICA component summaries.
- Save relevant outputs into `figure_source_data/` and/or `output_examples/` as needed.

### 7. IHC image preprocessing and measurement

- **Background harmonisation (optional but recommended for pale/low-stain tissue):**  
  - In Fiji, open `code/IHC_Background_Removal.ijm`.  
  - Inspect and, if needed, adjust the tunable parameters at the top of the script (hue/saturation/brightness thresholds, debug flags).  
  - Run the macro on a test image to confirm that the background is cleaned and the low-stain tissue is preserved.  
  - Use the batch mode (folder-level processing) if you wish to apply the same preprocessing to a set of images.

- **Batch cropping to standard ROIs:**  
  - Open `code/ImageJ_BatchCrop_SaveROIs.ijm` in Fiji.  
  - Specify input/output folders and confirm or update the fixed ROI coordinates in the dialog.  
  - Run the macro to generate cropped TIFF images with a consistent field-of-view.

- **Colour deconvolution and ROI measurements:**  
  - Open `code/ImageJ_ColorDeconvolution_BatchMeasurements.ijm`.  
  - Place the ROI sets (if applicable) in the expected location or load them into ROI Manager before running the macro.  
  - Run the macro on the folder containing cropped images; it will perform colour deconvolution, apply ROIs, and measure the chosen features.  
  - Export the ImageJ *Results* table as a tab-delimited file; this should match the structure documented in `IHC_ROI_measurements_header.tsv`.

### 8. Post-hoc BH and Fisher’s combined p-values (RNA-seq)

- If you wish to aggregate multiple tests across sites or related contrasts:  
  - Adapt `code/PostHoc_BH_and_Fisher_CombinedPvalues.R` to your result tables.  
  - Run the script to obtain BH-adjusted and Fisher-combined p-values per gene.

### 9. Gene annotation (RNA-seq)

- Whenever a results table contains Ensembl IDs but not gene symbols:  
  - Use `code/EnsemblID_to_GeneSymbol_Mapping_biomaRt.R` to add gene symbols via `biomaRt`.  
  - Export annotated tables to be used in final figures/tables and in `figure_source_data/`.

### 10. Generate figure panels (RNA-seq)

- After running the normalisation and differential expression workflows described above, use the figure-generation scripts in `code/` (see “Figure generation scripts (RNA-seq)”):  
  - `GTEx_Age_Lollipop_ByMatrisomeType.R` and `GTEx_Age_Lollipop_Combined_MatrisomeAndSites.R` to produce lollipop plots grouped by matrisome type and brain site (used in age-related regional panels).  
  - `Matrisome_Zscore_Heatmap_BrainSites.R` to generate Z-score heatmaps summarising regional matrisome expression patterns.  
  - `ROBO_Pathway_Age_ZscoreHeatmap.R` to generate Z-score heatmaps for the SLIT2–ROBO2 pathway across age groups.

All scripts contain in-script comments explaining the purpose of each major block of commands and the main design choices. Where the implementation details are complex, comments and section headers aim to make the logic explicit.

---

## Links to external methods repositories (if applicable)

Not applicable.

If accepted, this README will be updated to include author names, contact details, links to the preprint and final publication, and (where appropriate) links or references to full session information files for the main analyses.

---
