# Cell type deconvolution analysis workflow
- *requires imputed betas from differential methylation analysis*
1) Run cell-type deconvolution on DNAm data using EpiDISH to get cell type proportions: **cell_type_deconvolution_EpiDISH.R**
 
# Differential methylation analysis workflow
- *steps 4) and 5) require cell type proportions from cell_type_deconvolution_EpiDISH.R*
1) Impute betas from 450k array: **mliftover_impute_betas.R**
2) Run differential methylation analysis of DunedinPACE by genetic ancestry: **DML_DunedinPACE.R**
3) Run differential methylation analysis of Horvath clock by genetic ancestry: **DML_Horvath.R**
4) Run differential methylation analysis of DunedinPACE by genetic ancestry, covarying for cell type proportions: **DML_DunedinPACE_cell_type_props_adjusted.R**
5) Run differential methylation analysis of Horvath clock by genetic ancestry, covarying for cell type proportions: **DML_Horvath_cell_type_props_adjusted.R**
6) Analyse impact of cell type adjustment on differential methylation analyses: **DML_analyse_impact_of_cell_type_adjustment.R**

# Probe enrichment analysis workflow
 - *requires cpgs from differential methylation analysis*
1) **to be updated** Run gene enrichment analysis of differentially methylated cpgs using KnowYourCG and query enriched genes in GWAS Catalog: **gene_enrichment_analysis_with_GWASCatalog.R**
   - Use dev version of KnowYourCG and MSA platform
3) **to be updated** Run pathway analysis of genes enriched for differentially methylated cpgs: **pathway_analysis.R**
   - *requires genes from gene enrichment analysis*
   - Use dev version of KnowYourCG and MSA platform
4) **to be updated** Run probe enrichment analyses of differentially methylated probes: **probe_enrichment_analysis.R**
   - Use dev version of KnowYourCG and MSA platform
   - Add TI sig enrichment

# mQTL analysis workflow
1) Get genomic location of identified CpGs: **get_windows.R**
2) Locate imputed genotype data (chunked bed files): **locate_genotype_data_LPC.R**
3) Get ROIs from PLINK files for these data: **extract_ROI_from_PLINK_file.R**
4) Get population specific allele frequencies from gnomAD: **get_pop_specific_allele_freqs.R**
5) **to be updated** Run mQTL analysis: **run_mQTL_batch.R**
     - Update to fix error when using LD pruned and MAF filtered files
     - Update to include covariates in mQTL analysis
     - Color mQTL plots by ancestry group
6) **to be added** get summary table of SNPs from mQTL analysis with population-specific allele frequencies: **summarise_mQTL_allele_freqs.R**
7) **to be added** Sanity check for nearby common SNPs: 
