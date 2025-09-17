Differential methylation analysis workflow:
1) Impute betas from 450k array: mliftover_impute_betas.R
2) Run differentialy methylation analysis: DML_ancestry_test_sized_dot.R

Probe enrichment analysis workflow:
1) Run DMP by genetic ancestry: to be added
2) KnowYourCG analysis: to be added

mQTL analysis workflow:
1) Get genomic location of identified CpGs: get_windows.R
2) Locate imputed genotype data (chunked bed files): locate_genotype_data_LPC.R
3) Get ROIs from PLINK files for these data: extract_ROI_from_PLINK_file.R
4) Run mQTL analysis: to be added
     - EDIT SO THAT ASSOCIATIONS SIG AT FDR ARE CORRECTLY COUNTED AND PLOTTED
     - Color mQTL plots by ancestry group
     - run entirely within R (avoid command line steps)
     - secure password prompt
5) Sanity check for nearby common SNPs: to be added

Cell type deconvolution analysis workflow:
- requires imputed betas from differential methylation analysis workflow step #1
1) Get cell type proportions: cell_type_deconvolution_EpiDISH.R
2) Compare cell type proportions: compare_cell_type_proportions.R
