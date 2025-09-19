# 1) Differential methylation analysis workflow
a) Impute betas from 450k array: mliftover_impute_betas.R
b) Run differential methylation analysis: DML_ancestry_test_sized_dot.R

# 2) Probe enrichment analysis workflow
 - *requires cpgs from differential methylation analysis workflow step #1*
a) **to be added** Run probe enrichment on differentially methylated probes by genetic ancestry using KnowYourCG:

# 3) mQTL analysis workflow
a) Get genomic location of identified CpGs: get_windows.R
b) Locate imputed genotype data (chunked bed files): locate_genotype_data_LPC.R
c) Get ROIs from PLINK files for these data: extract_ROI_from_PLINK_file.R
d) **to be added** Get population specific allele frequencies from gnomAD: 
e) **to be updated** Run mQTL analysis: run_mQTL_batch.R
     - Update to fix error when using LD pruned and MAF filtered files
     - Update to include covariates in mQTL analysis
     - Update to include population-specific allele frequencies
     - Color mQTL plots by ancestry group
f) **to be added** Sanity check for nearby common SNPs: 

# 4) Cell type deconvolution analysis workflow
- *requires imputed betas from differential methylation analysis workflow step #1*
a) Get cell type proportions: cell_type_deconvolution_EpiDISH.R
b) Compare cell type proportions: compare_cell_type_proportions.R
