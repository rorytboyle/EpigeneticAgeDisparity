# This script locates PMBB genotype data (imputed genotypes as PLINK files) on the LPC
# Author: Rory Boyle & Nadia Dehghani rorytboyle@gmail.com
# Date: 17/09/2025

# chunk_files is output from locate_genotype_data_LPC.R
# results df includes:
# - CpG_ID: CpG identifier
# - CpG_Position: Position of CpG
# - Window_Distance: Window size used
# - variants_from_bim: Variants extracted from BIM file within the genomic window
# - variants_pruned: Number of variants removed by LD pruning
# - variants_removed_MAF: Number of variants removed by MAF filter
# - final_variants: Final number of variants in output


# Generate and execute PLINK commands with secure password prompt
run_plink_extractions <- function(chunk_files, 
                                  windows_df, 
                                  maf = 0.1, 
                                  ld_prune = TRUE,
                                  ld_window = 1500,
                                  ld_step = 50,
                                  ld_r2 = 0.8) {
  
  # Merge to get CpG positions
  merged <- merge(chunk_files, windows_df, by.x = "CpG_ID", by.y = "Input_CpG")
  
  # Initialize results tracking
  merged$variants_from_bim <- NA
  merged$variants_pruned <- NA
  merged$variants_removed_MAF <- NA
  merged$final_variants <- NA
  
  # Generate PLINK commands and parsing commands
  commands <- c()
  parse_commands <- c()
  
  for (i in 1:nrow(merged)) {
    cpg_id <- merged$CpG_ID[i]
    
    # Convert local mount path to LPC path
    local_path <- merged$File_Path[i]
    lpc_path <- gsub("~/lpc_mount/", "/static/PMBB/", local_path)
    lpc_path <- gsub("/Users/[^/]+/lpc_mount/", "/static/PMBB/", lpc_path)
    
    bfile <- gsub("\\.bed$", "", lpc_path)
    chromosome <- gsub("chr", "", merged$Chromosome[i])
    target_pos <- merged$CpG_Position[i]
    window_bp <- merged$Window_Distance[i]
    start_bp <- target_pos - window_bp
    end_bp <- target_pos + window_bp
    
    # Define output paths
    base_output <- paste0("/project/ftdc_external/PMBB/genetics/", cpg_id, "_", window_bp, "bp")
    
    if (ld_prune) {
      # Step 1: Extract region and perform LD pruning
      ld_output <- paste0(base_output, "_ld")
      cmd1 <- sprintf('echo "Processing %s - Step 1: LD pruning" && plink --bfile "%s" --chr %s --from-bp %d --to-bp %d --indep-pairwise %d %d %g --make-bed --out "%s" 2>&1 | tee "%s_step1.log"',
                      cpg_id, bfile, chromosome, start_bp, end_bp, ld_window, ld_step, ld_r2, ld_output, base_output)
      
      # Step 2: Extract independent SNPs and apply MAF filter
      final_output <- paste0(base_output, "_ld_maf")
      prune_file <- paste0(ld_output, ".prune.in")
      cmd2 <- sprintf('echo "Processing %s - Step 2: MAF filtering" && plink --bfile "%s" --extract "%s" --maf %g --recode A --out "%s" 2>&1 | tee "%s_step2.log"',
                      cpg_id, ld_output, prune_file, maf, final_output, base_output)
      
      commands <- c(commands, cmd1, cmd2)
      
      # Add parsing commands to extract statistics
      parse_cmd <- sprintf('
echo "=== Statistics for %s ===" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
echo "Step 1 (LD pruning):" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
grep -E "out of .* variants loaded|pass filters|Pruned .* variants|Pruning complete" "%s_step1.log" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
echo -n "Variants in prune.in: " >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
wc -l < "%s" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
echo "Step 2 (MAF filtering):" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
grep -E "variants loaded|--extract:|removed due to minor allele threshold|variants and.*people pass" "%s_step2.log" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
echo "---" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
', cpg_id, base_output, prune_file, base_output)
      
      parse_commands <- c(parse_commands, parse_cmd)
      
    } else {
      # If no LD pruning, just apply MAF filter directly
      output <- paste0(base_output, "_maf")
      cmd <- sprintf('echo "Processing %s - MAF filtering only" && plink --bfile "%s" --chr %s --from-bp %d --to-bp %d --maf %g --recode A --out "%s" 2>&1 | tee "%s.log"',
                     cpg_id, bfile, chromosome, start_bp, end_bp, maf, output, output)
      commands <- c(commands, cmd)
      
      # Add parsing command for non-LD version
      parse_cmd <- sprintf('
echo "=== Statistics for %s ===" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
grep -E "out of .* variants loaded|pass filters|removed due to minor allele threshold|variants and.*people pass" "%s.log" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
echo "---" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
', cpg_id, output)
      
      parse_commands <- c(parse_commands, parse_cmd)
    }
  }
  
  # Create script content
  script_content <- c(
    "#!/bin/bash",
    "module load plink/1.90Beta4.5",
    "",
    sprintf("# PLINK extraction script with MAF=%g and LD pruning=%s", maf, ifelse(ld_prune, "yes", "no")),
    if (ld_prune) sprintf("# LD parameters: window=%d kb, step=%d, r²=%g", ld_window, ld_step, ld_r2),
    "# Two-step process: 1) LD pruning, 2) Extract independent SNPs with MAF filter",
    "",
    "# Clear previous stats file",
    "rm -f /project/ftdc_external/PMBB/genetics/extraction_stats.txt",
    "",
    "# Run PLINK commands",
    commands,
    "",
    "# Parse statistics",
    parse_commands,
    "",
    "# Display summary",
    "echo '======= EXTRACTION SUMMARY ======='",
    "cat /project/ftdc_external/PMBB/genetics/extraction_stats.txt"
  )
  
  # Write script
  writeLines(script_content, "plink_extractions.sh")
  
  # Check if sshpass is available
  sshpass_check <- system("which sshpass", ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (sshpass_check != 0) {
    cat("sshpass not found. Please install it:\n")
    cat("Mac: brew install sshpass\n")
    cat("Linux: sudo apt-get install sshpass\n")
    cat("Or use manual approach: create_plink_script(chunk_files, windows_df)\n")
    return(NULL)
  }
  
  # Secure password prompt
  password <- get_secure_password()
  
  if (is.null(password)) {
    cat("Password entry cancelled or failed.\n")
    return(NULL)
  }
  
  # Transfer script
  cat("Transferring script...\n")
  scp_cmd <- sprintf('sshpass -p "%s" scp plink_extractions.sh rorytb@ftdcsub.pmacs.upenn.edu:~/', password)
  scp_result <- system(scp_cmd, ignore.stderr = FALSE)
  
  if (scp_result != 0) {
    cat("Failed to transfer script. Check password and connection.\n")
    # Clear password from memory before returning
    password <- NULL
    gc()
    return(NULL)
  }
  
  cat("Executing PLINK extractions on LPC...\n")
  ssh_cmd <- sprintf('sshpass -p "%s" ssh rorytb@ftdcsub.pmacs.upenn.edu "chmod +x plink_extractions.sh && bash plink_extractions.sh"', password)
  ssh_result <- system(ssh_cmd, ignore.stderr = FALSE)
  
  # Retrieve statistics file
  cat("Retrieving extraction statistics...\n")
  stats_file <- "/project/ftdc_external/PMBB/genetics/extraction_stats.txt"
  local_stats_file <- "extraction_stats.txt"
  scp_stats_cmd <- sprintf('sshpass -p "%s" scp rorytb@ftdcsub.pmacs.upenn.edu:%s %s', password, stats_file, local_stats_file)
  scp_stats_result <- system(scp_stats_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  # Clear password from memory after all commands are done
  password <- NULL
  gc()
  
  # Parse statistics if file was retrieved
  if (scp_stats_result == 0 && file.exists(local_stats_file)) {
    cat("\nParsing extraction statistics...\n")
    stats_content <- readLines(local_stats_file)
    
    # Parse stats for each CpG
    for (i in 1:nrow(merged)) {
      cpg_id <- merged$CpG_ID[i]
      cpg_stats <- extract_cpg_stats(stats_content, cpg_id, ld_prune = TRUE)
      
      merged$variants_from_bim[i] <- cpg_stats$variants_from_bim
      merged$variants_pruned[i] <- cpg_stats$variants_pruned
      merged$variants_removed_MAF[i] <- cpg_stats$variants_removed_MAF
      merged$final_variants[i] <- cpg_stats$final_variants
    }
    
    # Clean up local stats file
    unlink(local_stats_file)
    
    # Display summary
    cat("\n========= VARIANT STATISTICS SUMMARY =========\n")
    for (i in 1:nrow(merged)) {
      cat(sprintf("\n%s:\n", merged$CpG_ID[i]))
      cat(sprintf("  Variants from BIM file: %s\n", 
                  ifelse(is.na(merged$variants_from_bim[i]), "N/A", merged$variants_from_bim[i])))
      cat(sprintf("  Variants pruned (LD): %s\n", 
                  ifelse(is.na(merged$variants_pruned[i]), "N/A", merged$variants_pruned[i])))
      cat(sprintf("  Variants removed (MAF < %g): %s\n", maf,
                  ifelse(is.na(merged$variants_removed_MAF[i]), "N/A", merged$variants_removed_MAF[i])))
      cat(sprintf("  Final variants: %s\n", 
                  ifelse(is.na(merged$final_variants[i]), "N/A", merged$final_variants[i])))
    }
    cat("==============================================\n")
  } else {
    cat("Could not retrieve statistics file. Check logs on server for details.\n")
  }
  
  # Cleanup
  unlink("plink_extractions.sh")
  
  if (ssh_result == 0) {
    cat("\nPLINK extractions completed successfully!\n")
    cat(sprintf("Parameters used: MAF=%g, LD pruning=%s\n", maf, ifelse(ld_prune, "yes", "no")))
    if (ld_prune) {
      cat(sprintf("LD parameters: window=%d kb, step=%d, r²=%g\n", ld_window, ld_step, ld_r2))
    }
  } else {
    cat("Execution may have failed. Check output above.\n")
    cat("Note: You may need to run this in an interactive session (ibash) manually.\n")
  }
  
  # Return expanded results
  return_cols <- c("CpG_ID", "CpG_Position", "Window_Distance", 
                   "variants_from_bim", "variants_pruned", 
                   "variants_removed_MAF", "final_variants")
  return(merged[, return_cols])
}

# Function to parse statistics for a specific CpG
# helper: return first integer in a line (or NA)
extract_number <- function(line) {
  if (length(line) == 0) return(NA_real_)
  nums <- regmatches(line, gregexpr("\\d+", line))[[1]]
  if (length(nums) == 0) return(NA_real_)
  as.numeric(nums[1])
}

# robust parser for one CpG section of the stats file
extract_cpg_stats <- function(stats_content, cpg_id, ld_prune = TRUE) {
  stats <- list(
    variants_from_bim = NA_integer_,
    variants_pruned = NA_integer_,
    variants_removed_MAF = NA_integer_,
    final_variants = NA_integer_,
    after_pruning = NA_integer_  # helpful for debugging / checks
  )
  
  # locate section
  cpg_start <- grep(paste0("=== Statistics for ", cpg_id), stats_content)
  if (length(cpg_start) == 0) return(stats)
  cpg_end_candidates <- which(stats_content == "---")
  cpg_end <- cpg_end_candidates[cpg_end_candidates > cpg_start[1]][1]
  if (is.na(cpg_end)) cpg_end <- length(stats_content)
  cpg_section <- stats_content[cpg_start[1]:cpg_end]
  
  ## 1) variants_from_bim: first "variants loaded" line in the section
  loaded_lines <- grep("variants loaded", cpg_section, value = TRUE)
  if (length(loaded_lines) > 0) stats$variants_from_bim <- extract_number(loaded_lines[1])
  
  ## 2) after_pruning (prefer sources in this priority)
  prune_in <- grep("Variants in prune.in", cpg_section, value = TRUE)
  extract_line <- grep("--extract: .*variants remaining", cpg_section, value = TRUE)
  leaving_line <- grep("leaving \\d+", cpg_section, value = TRUE)
  pruned_line <- grep("^Pruned .* variants", cpg_section, value = TRUE)
  
  if (length(prune_in) > 0) {
    stats$after_pruning <- extract_number(prune_in[1])
  } else if (length(extract_line) > 0) {
    # often appears in step2 log (before MAF)
    stats$after_pruning <- extract_number(extract_line[1])
  } else if (length(leaving_line) > 0) {
    stats$after_pruning <- extract_number(leaving_line[1])
  } else if (length(pruned_line) > 0) {
    # try "leaving N" inside the Pruned line, otherwise infer from "Pruned X" and loaded
    leaving_m <- regmatches(pruned_line[1], regexpr("leaving \\d+", pruned_line[1]))
    if (length(leaving_m) && nzchar(leaving_m)) {
      stats$after_pruning <- as.numeric(sub("leaving ", "", leaving_m))
    } else {
      pruned_n <- extract_number(pruned_line[1])    # number pruned
      if (!is.na(stats$variants_from_bim) && !is.na(pruned_n)) {
        stats$after_pruning <- stats$variants_from_bim - pruned_n
        stats$variants_pruned <- pruned_n
      }
    }
  }
  
  ## 3) variants_pruned (if not already set)
  if (is.na(stats$variants_pruned) && !is.na(stats$variants_from_bim) && !is.na(stats$after_pruning)) {
    stats$variants_pruned <- stats$variants_from_bim - stats$after_pruning
  }
  
  ## 4) variants_removed_MAF: sum any matching lines (handles multiple lines if present)
  maf_lines <- grep("removed due to minor allele threshold", cpg_section, value = TRUE)
  if (length(maf_lines) > 0) {
    maf_nums <- sapply(maf_lines, extract_number)
    stats$variants_removed_MAF <- sum(maf_nums, na.rm = TRUE)
    if (is.na(stats$variants_removed_MAF)) stats$variants_removed_MAF <- NA_integer_
  }
  
  ## 5) final_variants: prefer the "X variants and Y people pass filters and QC" line (last occurrence)
  final_lines <- grep("variants and .*people pass filters and QC", cpg_section, value = TRUE)
  if (length(final_lines) > 0) {
    stats$final_variants <- extract_number(final_lines[length(final_lines)])
  } else if (!is.na(stats$after_pruning) && !is.na(stats$variants_removed_MAF)) {
    # fallback calculation if PLINK did not print the pass-line
    stats$final_variants <- stats$after_pruning - stats$variants_removed_MAF
  } else if (length(extract_line) > 0 && is.na(stats$final_variants)) {
    # last fallback: use --extract count (may be pre-MAF, so only a weak fallback)
    stats$final_variants <- extract_number(extract_line[length(extract_line)])
  }
  
  # sanity checks
  if (!is.na(stats$final_variants) && stats$final_variants < 0) stats$final_variants <- 0
  if (!is.na(stats$variants_pruned) && stats$variants_pruned < 0) stats$variants_pruned <- NA_integer_
  
  return(stats)
}


# Secure password prompt function
get_secure_password <- function() {
  # Try to use getPass package for secure password entry
  if (requireNamespace("getPass", quietly = TRUE)) {
    password <- getPass::getPass(msg = "Enter password for rorytb@ftdcsub.pmacs.upenn.edu: ")
  } else {
    cat("Note: Install 'getPass' package for more secure password entry:\n")
    cat("  install.packages('getPass')\n")
    cat("\nWARNING: Password may be visible in console history!\n")
    cat("Enter password for rorytb@ftdcsub.pmacs.upenn.edu: ")
    password <- readline()
  }
  
  return(password)
}

# Generate script for manual execution with ibash
create_plink_script <- function(chunk_files, 
                                windows_df, 
                                script_name = "plink_extractions.sh",
                                maf = 0.1,
                                ld_prune = TRUE,
                                ld_window = 1500,
                                ld_step = 50,
                                ld_r2 = 0.8,
                                return_stats = TRUE) {
  
  # Merge to get CpG positions
  merged <- merge(chunk_files, windows_df, by.x = "CpG_ID", by.y = "Input_CpG")
  
  # Initialize results tracking
  merged$variants_loaded <- NA
  merged$variants_after_region <- NA
  merged$variants_after_ld <- NA  # Number of independent SNPs after LD pruning
  merged$variants_removed_ld <- NA  # Number removed by LD pruning
  merged$variants_removed_maf <- NA
  merged$variants_final <- NA
  
  # Generate PLINK commands
  commands <- c()
  parse_commands <- c()
  
  for (i in 1:nrow(merged)) {
    cpg_id <- merged$CpG_ID[i]
    
    # Convert local mount path to LPC path
    local_path <- merged$File_Path[i]
    lpc_path <- gsub("~/lpc_mount/", "/static/", local_path)
    lpc_path <- gsub("/Users/[^/]+/lpc_mount/", "/static/", lpc_path)
    
    bfile <- gsub("\\.bed$", "", lpc_path)
    chromosome <- gsub("chr", "", merged$Chromosome[i])
    target_pos <- merged$CpG_Position[i]
    window_bp <- merged$Window_Distance[i]
    start_bp <- target_pos - window_bp
    end_bp <- target_pos + window_bp
    
    # Define output paths
    base_output <- paste0("/project/ftdc_external/PMBB/genetics/", cpg_id, "_", window_bp, "bp")
    
    if (ld_prune) {
      # Step 1: Extract region and perform LD pruning
      ld_output <- paste0(base_output, "_ld")
      cmd1 <- sprintf('echo "Processing %s - Step 1: LD pruning" && plink --bfile "%s" --chr %s --from-bp %d --to-bp %d --indep-pairwise %d %d %g --make-bed --out "%s" 2>&1 | tee "%s_step1.log"',
                      cpg_id, bfile, chromosome, start_bp, end_bp, ld_window, ld_step, ld_r2, ld_output, base_output)
      
      # Step 2: Extract independent SNPs and apply MAF filter
      final_output <- paste0(base_output, "_ld_maf")
      prune_file <- paste0(ld_output, ".prune.in")
      cmd2 <- sprintf('echo "Processing %s - Step 2: MAF filtering" && plink --bfile "%s" --extract "%s" --maf %g --recode A --out "%s" 2>&1 | tee "%s_step2.log"',
                      cpg_id, ld_output, prune_file, maf, final_output, base_output)
      
      commands <- c(commands, cmd1, cmd2)
      
      # Add parsing commands
      parse_cmd <- sprintf('
echo "=== Statistics for %s ===" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
echo "Step 1 (LD pruning):" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
grep -E "out of .* variants loaded|pass filters|Pruned .* variants|Pruning complete" "%s_step1.log" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
echo -n "Variants in prune.in: " >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
wc -l < "%s" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
echo "Step 2 (MAF filtering):" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
grep -E "variants loaded|--extract:|removed due to minor allele threshold|variants and.*people pass" "%s_step2.log" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
echo "---" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
', cpg_id, base_output, prune_file, base_output)
      
      parse_commands <- c(parse_commands, parse_cmd)
      
    } else {
      # If no LD pruning, just apply MAF filter directly
      output <- paste0(base_output, "_maf")
      cmd <- sprintf('echo "Processing %s - MAF filtering only" && plink --bfile "%s" --chr %s --from-bp %d --to-bp %d --maf %g --recode A --out "%s" 2>&1 | tee "%s.log"',
                     cpg_id, bfile, chromosome, start_bp, end_bp, maf, output, output)
      commands <- c(commands, cmd)
      
      # Add parsing command
      parse_cmd <- sprintf('
echo "=== Statistics for %s ===" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
grep -E "out of .* variants loaded|pass filters|removed due to minor allele threshold|variants and.*people pass" "%s.log" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
echo "---" >> /project/ftdc_external/PMBB/genetics/extraction_stats.txt
', cpg_id, output)
      
      parse_commands <- c(parse_commands, parse_cmd)
    }
  }
  
  # Create script content
  script_content <- c(
    "#!/bin/bash",
    "module load plink/1.90Beta4.5",
    "",
    sprintf("# PLINK extraction script with MAF=%g and LD pruning=%s", maf, ifelse(ld_prune, "yes", "no")),
    if (ld_prune) sprintf("# LD parameters: window=%d kb, step=%d, r²=%g", ld_window, ld_step, ld_r2),
    if (ld_prune) "# Two-step process: 1) LD pruning, 2) Extract independent SNPs with MAF filter",
    "",
    "# Clear previous stats file",
    "rm -f /project/ftdc_external/PMBB/genetics/extraction_stats.txt",
    "",
    "# Run PLINK commands",
    commands,
    "",
    "# Parse statistics",
    parse_commands,
    "",
    "# Display summary",
    "echo '======= EXTRACTION SUMMARY ======='",
    "cat /project/ftdc_external/PMBB/genetics/extraction_stats.txt"
  )
  
  # Write script
  writeLines(script_content, script_name)
  system(paste("chmod +x", script_name))
  
  cat(paste("Created script:", script_name, "\n"))
  cat(sprintf("Parameters: MAF=%g, LD pruning=%s\n", maf, ifelse(ld_prune, "yes", "no")))
  if (ld_prune) {
    cat(sprintf("LD parameters: window=%d kb (%.1f Mb), step=%d, r²=%g\n", 
                ld_window, ld_window/1000, ld_step, ld_r2))
    cat("Two-step process:\n")
    cat("  Step 1: Extract region and perform LD pruning (creates .prune.in file)\n")
    cat("  Step 2: Extract independent SNPs from .prune.in and apply MAF filter\n")
  }
  cat("\nTo run on LPC:\n")
  cat(paste("1. scp", script_name, "rorytb@ftdcsub.pmacs.upenn.edu:~/\n"))
  cat("2. ssh rorytb@ftdcsub.pmacs.upenn.edu\n")
  cat("3. ibash\n")
  cat(paste("4. bash", script_name, "\n"))
  
  if (return_stats) {
    cat("\nTo retrieve statistics after running:\n")
    cat("scp rorytb@ftdcsub.pmacs.upenn.edu:/project/ftdc_external/PMBB/genetics/extraction_stats.txt .\n")
  }
  
  # Return expanded results columns
  return_cols <- c("CpG_ID", "CpG_Position", "Window_Distance", 
                   "variants_from_bim", "variants_pruned", 
                   "variants_removed_MAF", "final_variants")
  return(merged[, return_cols])
}

# Usage ####
results <- run_plink_extractions(chunk_files, windows, maf = 0.1, ld_prune = TRUE, ld_window = 1500, ld_step = 50, ld_r2 = 0.8)

# Example usage with custom parameters
# results <- run_plink_extractions(chunk_files, windows, 
#                                 maf = 0.15, 
#                                 ld_r2 = 0.9,
#                                 ld_window = 2000)

# Example without LD pruning (single-step process with just MAF filter)
# results <- run_plink_extractions(chunk_files, windows, 
#                                 maf = 0.1, 
#                                 ld_prune = FALSE)