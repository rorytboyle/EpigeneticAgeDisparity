library(DunedinPACE)

cpg_sites = getRequiredProbes()
cpg_sites=cpg_sites$DunedinPACE
saveRDS(cpg_sites, file = "20250212_cpg_DunedinPACE.rds")

cpg_sites =  getRequiredProbes(backgroundList = TRUE)
cpg_sites=cpg_sites$DunedinPACE
saveRDS(cpg_sites, file = "20250212_cpg_DunedinPACE_all.rds")