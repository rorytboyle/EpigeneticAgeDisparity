# This script runs mQTL analyses for provided CpGs using MatrixEQTL
# Author: Nadia Dehghani & Rory Boyle rorytboyle@gmail.com
# Date: 10/09/2025
# https://github.com/andreyshabalin/MatrixEQTL

## libraries --------
library(tidyverse)
library(dplyr)
library(MatrixEQTL)
library(ggplot2)
library(gridExtra)
library(getPass)  

## Remote connection functions --------
mount_remote_directory <- function(remote_user = "rorytb", 
                                   remote_host = "ftdcsub1",
                                   remote_path = "/project/ftdc_external/PMBB/",
                                   local_mount = "/Users/rorytb/PMACS_remote") {
  
  cat("Setting up remote connection...\n")
  
  # Create local mount directory if it doesn't exist
  if (!dir.exists(local_mount)) {
    dir.create(local_mount, recursive = TRUE)
  }
  
  # Check if already mounted
  mount_check <- system(paste("mount | grep", local_mount), intern = TRUE)
  if (length(mount_check) > 0) {
    cat("Directory already mounted.\n")
    return(TRUE)
  }
  
  # Get password securely
  password <- getPass::getPass(msg = paste("Enter password for", remote_user, "@", remote_host, ": "))
  
  # Create temporary expect script for password automation
  expect_script <- tempfile(fileext = ".exp")
  cat(sprintf('#!/usr/bin/expect
set timeout 20
spawn sshfs %s@%s:%s %s -o defer_permissions,volname=project
expect "password:"
send "%s\\r"
expect eof
', remote_user, remote_host, remote_path, local_mount, password), file = expect_script)
  
  # Make script executable and run it
  system(paste("chmod +x", expect_script))
  mount_result <- system(paste("expect", expect_script))
  
  # Clean up temporary script
  unlink(expect_script)
  
  # Clear password from memory
  rm(password)
  gc()
  
  if (mount_result == 0) {
    cat("Successfully mounted remote directory.\n")
    return(TRUE)
  } else {
    cat("Failed to mount remote directory.\n")
    return(FALSE)
  }
}

unmount_remote_directory <- function(local_mount = "/Users/rorytb/PMACS_remote") {
  cat("Unmounting remote directory...\n")
  
  # Try different unmount commands depending on OS
  if (Sys.info()["sysname"] == "Darwin") {
    # macOS
    unmount_result <- system(paste("diskutil umount force", local_mount))
  } else {
    # Linux/Unix
    unmount_result <- system(paste("fusermount -u", local_mount))
  }
  
  if (unmount_result == 0) {
    cat("Successfully unmounted remote directory.\n")
    return(TRUE)
  } else {
    cat("Failed to unmount remote directory. You may need to unmount manually.\n")
    return(FALSE)
  }
}

## Main mQTL analysis function --------
run_mqtl_analysis <- function(cpg_info, betas_data, base_output_dir, base_genetics_dir) {
  
  cpg_id <- cpg_info$Input_CpG
  cpg_position <- cpg_info$CpG_Position
  window_distance <- cpg_info$Window_Distance
  
  cat("\n" , rep("=", 60), "\n")
  cat("Processing CpG:", cpg_id, "\n")
  cat("Position:", cpg_info$Chromosome, ":", cpg_position, "\n")
  cat("Window: +/-", window_distance, "bp\n")
  cat(rep("=", 60), "\n")
  
  # Setup paths
  plink_file <- file.path(base_genetics_dir, paste0(cpg_id, "_", window_distance, "bp_ld_maf.raw"))
  output_dir <- file.path(base_output_dir, paste0("mQTL_", cpg_id))
  
  # Check if PLINK file exists
  if (!file.exists(plink_file)) {
    cat("ERROR: PLINK file not found:", plink_file, "\n")
    return(list(cpg_id = cpg_id, status = "failed", error = "PLINK file missing"))
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  tryCatch({
    
    ## Load and process genotype data --------
    cat("Loading genotype data...\n")
    geno <- read.table(plink_file, header = TRUE, sep = " ")
    geno$IID <- as.character(geno$IID)
    geno <- subset(geno, select = -c(FID, PAT, MAT, SEX, PHENOTYPE))
    
    if (ncol(geno) <= 1) {
      cat("ERROR: No SNPs found in PLINK file\n")
      return(list(cpg_id = cpg_id, status = "failed", error = "No SNPs"))
    }
    
    cat("Found", ncol(geno) - 1, "SNPs\n")
    
    ## Process methylation data --------
    cat("Processing methylation data...\n")
    
    # Find CpG in methylation data (handle suffix matching)
    cpg_manifest_id <- cpg_info$Manifest_ID
    if (!cpg_manifest_id %in% rownames(betas_data)) {
      cat("ERROR: CpG", cpg_manifest_id, "not found in methylation data\n")
      return(list(cpg_id = cpg_id, status = "failed", error = "CpG not in methylation data"))
    }
    
    # Extract methylation data for this CpG
    betas_subset <- betas_data[cpg_manifest_id, , drop = FALSE]
    betas2 <- rownames_to_column(as.data.frame(t(betas_subset)), var = "IID")
    names(betas2)[2] <- cpg_manifest_id
    
    ## Merge genotype and methylation data --------
    cat("Merging datasets...\n")
    INDD_betas_geno <- merge(geno, betas2, by = "IID")
    
    if (nrow(INDD_betas_geno) == 0) {
      cat("ERROR: No overlapping samples between genotype and methylation data\n")
      return(list(cpg_id = cpg_id, status = "failed", error = "No overlapping samples"))
    }
    
    cat("Found", nrow(INDD_betas_geno), "overlapping samples\n")
    
    ## Prepare data for MatrixEQTL --------
    # Genotype matrix
    INDD_geno <- unique(INDD_betas_geno[, 1:(ncol(INDD_betas_geno) - 1)])
    n <- INDD_geno$IID
    tgeno <- as.data.frame(t(INDD_geno[, -1]))
    colnames(tgeno) <- n
    SNPs <- rownames_to_column(tgeno, "snpid")
    SNPs <- na.omit(SNPs)
    
    # Methylation matrix
    INDD_betas <- unique(subset(INDD_betas_geno, select = c(1, ncol(INDD_betas_geno))))
    n <- INDD_betas$IID
    tINDD_beta_CG <- as.data.frame(t(INDD_betas[, -1]))
    colnames(tINDD_beta_CG) <- n
    GE <- rownames_to_column(tINDD_beta_CG, "geneid")
    GE$geneid <- cpg_manifest_id
    
    # Empty covariates
    Covariates <- data.frame(id = character(0))
    
    ## Export data files --------
    write.table(SNPs, file.path(output_dir, "SNPs.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(GE, file.path(output_dir, "GE.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(Covariates, file.path(output_dir, "Covariates.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    
    ## Run MatrixEQTL --------
    cat("Running MatrixEQTL analysis...\n")
    
    useModel = modelLINEAR
    output_file_name = file.path(output_dir, "PMBB_mQTL_output")
    pvOutputThreshold = 1
    errorCovariance = numeric()
    
    # Load data
    snps = SlicedData$new()
    snps$fileDelimiter = "\t"
    snps$fileOmitCharacters = "NA"
    snps$fileSkipRows = 1
    snps$fileSkipColumns = 1
    snps$fileSliceSize = 2000
    snps$LoadFile(file.path(output_dir, "SNPs.txt"))
    
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"
    gene$fileOmitCharacters = "NA"
    gene$fileSkipRows = 1
    gene$fileSkipColumns = 1
    gene$fileSliceSize = 2000
    gene$LoadFile(file.path(output_dir, "GE.txt"))
    
    if (nrow(SNPs) == 0) {
      cat("No SNPs remaining after filtering. Skipping CpG", cpg_id, "\n")
      return(list(cpg_id = cpg_id, status = "failed", error = "No SNPs after filtering"))
    }
    
    if (nrow(GE) == 0) {
      cat("No methylation data for CpG", cpg_id, "\n")
      return(list(cpg_id = cpg_id, status = "failed", error = "No methylation data"))
    }
    
    # Run analysis
    me = Matrix_eQTL_engine(
      snps = snps,
      gene = gene,
      output_file_name = output_file_name,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = FALSE,  # Reduce output for batch processing
      pvalue.hist = FALSE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE
    )
    
    ## Process results --------
    mqtl_results <- me$all$eqtls
    
    if (nrow(mqtl_results) == 0) {
      cat("WARNING: No associations found\n")
      return(list(cpg_id = cpg_id, status = "completed", n_snps = 0, n_significant = 0))
    }
    
    # Handle FDR column naming (MatrixEQTL uses uppercase 'FDR')
    if ("FDR" %in% names(mqtl_results) && !"fdr" %in% names(mqtl_results)) {
      mqtl_results$fdr <- mqtl_results$FDR
    }
    
    # Calculate genomic distances
    calculate_genomic_distance <- function(snp_names, cpg_pos) {
      snp_positions <- as.numeric(gsub(".*_([0-9]+)_[A-Z]+_[A-Z]+.*", "\\1", snp_names))
      distances <- abs(snp_positions - cpg_pos)
      return(distances)
    }
    
    mqtl_results$genomic_distance <- calculate_genomic_distance(mqtl_results$snps, cpg_position)
    mqtl_results <- mqtl_results[order(mqtl_results$pvalue), ]
    
    # Add CpG info to results
    mqtl_results$cpg_id <- cpg_id
    mqtl_results$cpg_position <- cpg_position
    
    ## Generate plots for significant associations --------
    significance_threshold <- 0.05
    significant_hits <- mqtl_results[mqtl_results$pvalue < significance_threshold, ]
    
    cat("Found", nrow(significant_hits), "significant associations (P <", significance_threshold, ")\n")
    
    # Create distance vs association strength plot
    cat("Creating distance vs association plot...\n")
    
    # Define significance thresholds
    fdr_threshold <- 0.05
    nominal_threshold <- 0.05
    
    # Create significance categories for cleaner plotting
    mqtl_results$significance <- "Non-significant"
    mqtl_results$significance[mqtl_results$pvalue < nominal_threshold] <- "Nominal"
    if ("fdr" %in% names(mqtl_results)) {
      mqtl_results$significance[mqtl_results$fdr < fdr_threshold] <- "FDR"
    }
    mqtl_results$significance <- factor(mqtl_results$significance,
                                        levels = c("Non-significant", "Nominal", "FDR"))
    
    # Count significant associations
    n_fdr_significant <- 0
    if ("fdr" %in% names(mqtl_results)) {
      n_fdr_significant <- sum(mqtl_results$fdr < fdr_threshold, na.rm = TRUE)
    }
    n_nominal_significant <- sum(mqtl_results$pvalue < nominal_threshold, na.rm = TRUE)
    
    cat("FDR significant SNPs:", n_fdr_significant, "\n")
    if (n_fdr_significant > 0) {
      fdr_values <- mqtl_results$fdr[mqtl_results$fdr < fdr_threshold]
      cat("FDR values range:", min(fdr_values, na.rm = TRUE), "to", max(fdr_values, na.rm = TRUE), "\n")
    }
    
    distance_plot <- ggplot(mqtl_results, aes(x = genomic_distance, y = -log10(pvalue))) +
      geom_hline(yintercept = -log10(nominal_threshold), linetype = "dashed", color = "red", alpha = 0.7) +
      geom_point(aes(color = significance, fill = significance, size = significance, shape = significance), alpha = 0.8) +
      scale_color_manual(values = c("Non-significant" = "gray60", "Nominal" = "red", "FDR" = "darkred"),
                         name = "Significance") +
      scale_fill_manual(values = c("Non-significant" = "gray60", "Nominal" = "red", "FDR" = "yellow"),
                        name = "Significance") +
      scale_size_manual(values = c("Non-significant" = 1.5, "Nominal" = 1.5, "FDR" = 2.5),
                        name = "Significance") +
      scale_shape_manual(values = c("Non-significant" = 16, "Nominal" = 16, "FDR" = 21),
                         name = "Significance") +
      theme_minimal() +
      theme(legend.position = "bottom") +
      xlab("Distance from CpG (bp)") +
      ylab("-log10(P-value)") +
      ggtitle(paste("mQTL Association Strength vs Distance:", cpg_id),
              subtitle = paste("n =", nrow(mqtl_results), "SNPs,", n_nominal_significant, "nominal sig,", n_fdr_significant, "FDR sig")) +
      scale_x_continuous(labels = scales::comma) +
      annotate("text", x = Inf, y = -log10(nominal_threshold),
               label = paste("P =", nominal_threshold),
               hjust = 1.1, vjust = -0.5, size = 3, color = "red")
    
    # Add annotation for FDR significant count
    if (n_fdr_significant > 0) {
      distance_plot <- distance_plot +
        annotate("text", x = -Inf, y = Inf,
                 label = paste("FDR < 0.05:", n_fdr_significant, "SNPs"),
                 hjust = -0.1, vjust = 1.5, size = 3, color = "darkred", fontface = "bold")
    }
    
    # Save distance plot
    distance_filename <- file.path(output_dir, "distance_vs_association.png")
    ggsave(distance_filename, plot = distance_plot, width = 10, height = 6, dpi = 300)
    
    # Display distance plot in output
    print(distance_plot)
    
    if (nrow(significant_hits) > 0) {
      
      # Create plots for significant hits
      plots_list <- list()
      
      for (i in 1:min(nrow(significant_hits), 9)) {  # Limit to top 9
        hit <- significant_hits[i, ]
        snp_name <- hit$snps
        
        if (snp_name %in% names(INDD_betas_geno)) {
          
          # Create plot
          plot_data <- INDD_betas_geno
          plot_data[[snp_name]] <- as.factor(plot_data[[snp_name]])
          plot_data[[cpg_manifest_id]] <- as.numeric(plot_data[[cpg_manifest_id]])
          
          p <- ggplot(plot_data, aes_string(x = snp_name, y = cpg_manifest_id)) + 
            geom_violin(width = 1, alpha = 0.3) + 
            geom_jitter(position = position_jitter(width = 0.2), size = 0.5) + 
            geom_boxplot(color = 'seagreen3', coef = Inf, width = 0.3, alpha = 0.7) + 
            theme_minimal() +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                  axis.title = element_text(size = 10),
                  plot.title = element_text(size = 11)) + 
            xlab(paste("Genotype:", snp_name)) +
            ylab(paste("Methylation:", cpg_id)) +
            ggtitle(paste0("mQTL #", i, ": ", cpg_id),
                    subtitle = paste0("P = ", format(hit$pvalue, scientific = TRUE, digits = 2),
                                      ", Dist = ", format(hit$genomic_distance, big.mark = ","), " bp"))
          
          plots_list[[i]] <- p
          
          # Save individual plot
          plot_filename <- file.path(output_dir, paste0("mQTL_plot_", i, "_",
                                                        gsub("[^A-Za-z0-9_]", "_", snp_name), ".png"))
          ggsave(plot_filename, plot = p, width = 8, height = 6, dpi = 300)
        }
      }
      
      # Create combined plot
      if (length(plots_list) > 1) {
        n_plots <- length(plots_list)
        n_cols <- min(3, ceiling(sqrt(n_plots)))
        n_rows <- ceiling(n_plots / n_cols)
        
        combined_plot <- do.call(grid.arrange, c(plots_list, ncol = n_cols, nrow = n_rows))
        
        combined_filename <- file.path(output_dir, paste0("mQTL_combined_plots.png"))
        ggsave(combined_filename, plot = combined_plot,
               width = 4 * n_cols, height = 4 * n_rows, dpi = 300)
      }
    }
    
    # Save results
    results_filename <- file.path(output_dir, "mqtl_results_with_distances.csv")
    write.csv(mqtl_results, results_filename, row.names = FALSE)
    
    # Create summary statistics
    summary_stats <- list(
      cpg_id = cpg_id,
      status = "completed",
      n_snps = nrow(mqtl_results),
      n_significant = nrow(significant_hits),
      n_fdr_significant = n_fdr_significant,
      n_samples = nrow(INDD_betas_geno),
      top_pvalue = min(mqtl_results$pvalue),
      min_distance = min(mqtl_results$genomic_distance),
      max_distance = max(mqtl_results$genomic_distance)
    )
    
    if (nrow(significant_hits) > 0) {
      top_hit <- significant_hits[1, ]
      summary_stats$top_snp <- top_hit$snps
      summary_stats$top_beta <- top_hit$beta
      summary_stats$top_distance <- top_hit$genomic_distance
    }
    
    cat("Analysis completed successfully!\n")
    return(summary_stats)
    
  }, error = function(e) {
    cat("ERROR during analysis:", e$message, "\n")
    return(list(cpg_id = cpg_id, status = "failed", error = e$message))
  })
}

## Main execution function --------
run_batch_mqtl_analysis <- function(windows_data,
                                    betas_file_path = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/20250710_raw_betas_811_DunedinPACE_bw_disparity.rds",
                                    base_genetics_dir = "/Users/rorytb/PMACS_remote/genetics",
                                    base_output_dir = "/Users/rorytb/PMACS_remote/genetics",
                                    use_remote = TRUE,
                                    remote_user = "rorytb",
                                    remote_host = "ftdcsub1",
                                    remote_path = "/project/ftdc_external/PMBB/",
                                    local_mount = "/Users/rorytb/PMACS_remote") {
  
  cat("Starting batch mQTL analysis for", nrow(windows_data), "CpGs...\n")
  
  # Mount remote directory if needed
  if (use_remote) {
    mount_success <- mount_remote_directory(
      remote_user = remote_user,
      remote_host = remote_host,
      remote_path = remote_path,
      local_mount = local_mount
    )
    
    if (!mount_success) {
      stop("Failed to mount remote directory. Please check your connection settings.")
    }
  }
  
  cat("Loading methylation data...\n")
  
  # Load methylation data once
  betas <- readRDS(betas_file_path)
  
  # Initialize results storage
  all_results <- list()
  summary_table <- data.frame()
  
  # Process each CpG
  for (i in 1:nrow(windows_data)) {
    
    cpg_info <- windows_data[i, ]
    
    result <- run_mqtl_analysis(
      cpg_info = cpg_info,
      betas_data = betas,
      base_output_dir = base_output_dir,
      base_genetics_dir = base_genetics_dir
    )
    
    all_results[[i]] <- result
    if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
      summary_table <- bind_rows(summary_table, as.data.frame(result))
    }    
    # Progress update
    cat("Completed", i, "of", nrow(windows_data), "CpGs\n")
  }
  
  ## Create overall summary --------
  summary_filename <- file.path(base_output_dir, "batch_mqtl_summary.csv")
  write.csv(summary_table, summary_filename, row.names = FALSE)
  
  cat("\n" , rep("=", 60), "\n")
  cat("BATCH ANALYSIS COMPLETE\n")
  cat(rep("=", 60), "\n")
  cat("Total CpGs processed:", nrow(summary_table), "\n")
  cat("Successful analyses:", sum(summary_table$status == "completed"), "\n")
  cat("Failed analyses:", sum(summary_table$status == "failed"), "\n")
  
  if (sum(summary_table$status == "completed") > 0) {
    completed <- summary_table[summary_table$status == "completed", ]
    cat("Total significant associations:", sum(completed$n_significant, na.rm = TRUE), "\n")
    cat("Total FDR-significant associations:", sum(completed$n_fdr_significant, na.rm = TRUE), "\n")
    cat("CpGs with significant hits:", sum(completed$n_significant > 0, na.rm = TRUE), "\n")
    cat("CpGs with FDR-significant hits:", sum(completed$n_fdr_significant > 0, na.rm = TRUE), "\n")
    cat("Summary saved to:", summary_filename, "\n")
  }
  
  # Unmount remote directory if it was mounted
  if (use_remote) {
    unmount_success <- unmount_remote_directory(local_mount)
    if (!unmount_success) {
      cat("Note: You may need to manually unmount the remote directory.\n")
    }
  }
  
  return(list(results = all_results, summary = summary_table))
}

## Usage --------
# Run analysis with secure password prompt
batch_results <- run_batch_mqtl_analysis(
  windows_data = windows,
  betas_file_path = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/20250710_raw_betas_811_DunedinPACE_bw_disparity.rds",
  base_genetics_dir = "/Users/rorytb/PMACS_remote/genetics",
  base_output_dir = "/Users/rorytb/PMACS_remote/genetics",
  use_remote = TRUE,
  remote_user = "rorytb",
  remote_host = "ftdcsub1",
  remote_path = "/project/ftdc_external/PMBB/",
  local_mount = "/Users/rorytb/PMACS_remote"
)

# Or if already mounted/using local files, set use_remote = FALSE
# batch_results <- run_batch_mqtl_analysis(
#   windows_data = windows,
#   betas_file_path = "/path/to/betas.rds",
#   base_genetics_dir = "/path/to/genetics",
#   base_output_dir = "/path/to/output",
#   use_remote = FALSE
# )
