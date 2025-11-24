betas <- readRDS("20250529_raw_betas_848.rds")
rownames(betas) <- sub("_.*", "", rownames(betas))

horvath_df=readRDS("/home/xuh6/20240920_Hao/projects/PMBB_data_analysis/20251124_cpg_Horvath_weights.rds")
common_probes <- intersect(rownames(betas), horvath_df$CpG)
length(common_probes)  #  347
head(common_probes)  # Check the first few common probes
# save common_probes to a rds file
saveRDS(common_probes, "20251124_cpg_Horvath_intersect_MSA.rds")




