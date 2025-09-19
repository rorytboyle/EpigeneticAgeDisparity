# This script runs probe enrichment analyses for CpGs from a differential methylation analysis
# Author: Rory Boyle & Nadia Dehghani rorytboyle@gmail.com
# Date: 19/09/2025

# Cache available databases - could reduce load by specifying priority databases below
# sesameData::sesameDataCache(data_titles=
# c("EPIC.address","genomeInfo.hg38","probeIDSignature"))
sesameDataCache()

# Read in CpG list from differential methylation results
dmp_cgs <- read.csv('/Users/rorytb/Library/CloudStorage/Box-Box/PMBB for Rory/20250819_171_stats.csv')

# Get CpG list from CpGs that are hypermethylated in African ancestry beyond FDR sig
dmp_cpgs_hyper_AA <- dmp_cgs %>%
  filter(adj.P.Val < 0.05) %>%
  filter(Ancestry == "Higher in AFR")
query_cpgs_hyper_AA <- dmp_cpgs_hyper_AA$CpG

# Check what's available
available_dbs <- listDBGroups("EPIC")

# Basic enrichment analyses ####
# Get informative databases for hypermethylation analysis
priority_databases <- c(
  "KYCG.EPIC.chromHMM.20211020",        # Chromatin states
  "KYCG.EPIC.TFBSconsensus.20211013",            # Transcription factor binding sites  
  "KYCG.EPIC.CGI.20210713",              # CpG islands
  "KYCG.EPIC.HMconsensus.20211013"
)

# chrom_res <- testEnrichment(query_cpgs_hyper_AA, platform="EPIC", databases = "KYCG.EPIC.chromHMM.20211020")
tfbs_res <- testEnrichment(query_cpgs_hyper_AA, platform="EPIC", databases = "KYCG.EPIC.TFBSconsensus.20211013") 
cgi_res <- testEnrichment(query_cpgs_hyper_AA, platform="EPIC", databases = "KYCG.EPIC.CGI.20210713")
# hist_res <- testEnrichment(query_cpgs_hyper_AA, platform="EPIC", databases = "KYCG.EPIC.HMconsensus.20211013")

# Gene enrichment analyses ####
# get dbs
dbs <- buildGeneDBs(query_cpgs_hyper_AA, platform = "EPIC")
gene_res <- testEnrichment(query_cpgs_hyper_AA, dbs, platform="EPIC")

# Visualize and interpret results TO BE ADDED ####

# Pathway analysis ####
library(gprofiler2)

sig_genes <- gene_res$gene_name[gene_res$FDR < 0.05]
pathway_res <- gost(sig_genes, organism = "hsapiens")
pathway_res$result
gostplot(pathway_res, capped = FALSE, interactive = FALSE)
