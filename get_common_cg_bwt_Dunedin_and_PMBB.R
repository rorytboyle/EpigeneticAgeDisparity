betas <- readRDS("20250529_raw_betas_848.rds")

rownames(betas) <- sub("_.*", "", rownames(betas))
pace=readRDS("/home/xuh6/20240920_Hao/projects/PMBB_data_analysis/20250212_cpg_DunedinPACE.rds")
pace_all=readRDS("/home/xuh6/20240920_Hao/projects/PMBB_data_analysis/20250212_cpg_DunedinPACE_all.rds")
common_probes <- intersect(rownames(betas), pace)
length(common_probes)  # Number of common probes
common_probes_all <- intersect(rownames(betas), pace_all)
length(common_probes_all)  # Number of common probes in the all set
# save common_probes to a rds file
saveRDS(common_probes, "20250529_cpg_DunedinPACE_intersect_MSA.rds")
saveRDS(common_probes_all, "20250529_cpg_DunedinPACE_all_intersect_MSA.rds")
common_probes <- readRDS("20250529_cpg_DunedinPACE_intersect_MSA.rds")
length(common_probes)  # Number of common probes
head(common_probes)  # Check the first few common probes
common_probes_all <- readRDS("20250529_cpg_DunedinPACE_all_intersect_MSA.rds")
length(common_probes_all)  # Number of common probes in the all set
head(common_probes_all)  # Check the first few common probes in the all set
