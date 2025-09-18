# Differential methylation analysis workflow
1) Impute betas from 450k array: mliftover_impute_betas.R
2) Run differential methylation analysis: DML_ancestry_test_sized_dot.R

# Probe enrichment analysis workflow
1) **to be added** Run DMP by genetic ancestry:
2) **to be added** KnowYourCG analysis:

# mQTL analysis workflow
1) Get genomic location of identified CpGs: get_windows.R
2) Locate imputed genotype data (chunked bed files): locate_genotype_data_LPC.R
3) Get ROIs from PLINK files for these data: extract_ROI_from_PLINK_file.R
4) **to be updated** Run mQTL analysis: run_mQTL_batch.R
     - Update to fix error when using LD pruned and MAF filtered files
     - Update to include covariates in mQTL analysis
     - Color mQTL plots by ancestry group
5) **to be added** Sanity check for nearby common SNPs: 

# Cell type deconvolution analysis workflow
- *requires imputed betas from differential methylation analysis workflow step #1*
1) Get cell type proportions: cell_type_deconvolution_EpiDISH.R
2) Compare cell type proportions: compare_cell_type_proportions.R
