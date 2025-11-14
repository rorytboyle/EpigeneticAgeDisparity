# Differential methylation analysis workflow
1) Impute betas from 450k array: mliftover_impute_betas.R
2) Run differential methylation analysis by genetic ancestry: DML_ancestry_test_sized_dot.R

# Probe enrichment analysis workflow
 - *requires cpgs from differential methylation analysis*
1) Run gene enrichment analysis of differentially methylated probes using KnowYourCG and query enriched genes in GWAS Catalog: gene_enrichment_analysis_with_GWASCatalog.R
2) **to be added** Run location analysis of differentially methylated probes: genomic_location_analysis.R
   - see probe_enrichment_analysis_WORKING.R FOR OTHER STOPS IN ENRICHMENT ANALYSIS

# mQTL analysis workflow
1) Get genomic location of identified CpGs: get_windows.R
2) Locate imputed genotype data (chunked bed files): locate_genotype_data_LPC.R
3) Get ROIs from PLINK files for these data: extract_ROI_from_PLINK_file.R
4) Get population specific allele frequencies from gnomAD: get_pop_specific_allele_freqs.R
5) **to be updated** Run mQTL analysis: run_mQTL_batch.R
     - Update to fix error when using LD pruned and MAF filtered files
     - Update to include covariates in mQTL analysis
     - Color mQTL plots by ancestry group
6) **to be added** get summary table of SNPs from mQTL analysis with population-specific allele frequencies: summarise_mQTL_allele_freqs.R
7) **to be added** Sanity check for nearby common SNPs: 

# Cell type deconvolution analysis workflow
- *requires imputed betas from differential methylation analysis*
1) Get cell type proportions: cell_type_deconvolution_EpiDISH.R
2) Compare cell type proportions: compare_cell_type_proportions.R
