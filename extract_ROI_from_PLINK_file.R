# Generate and execute PLINK commands with password prompt
run_plink_extractions <- function(chunk_files, windows_df) {
  
  # Merge to get CpG positions
  merged <- merge(chunk_files, windows_df, by.x = "CpG_ID", by.y = "Input_CpG")
  
  # Generate PLINK commands
  commands <- c()
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
    output <- paste0("/project/ftdc_external/PMBB/genetics/", cpg_id, "_", window_bp, "bp")
    
    cmd <- sprintf('plink --bfile "%s" --chr %s --from-bp %d --to-bp %d --recode A --out "%s"',
                   bfile, chromosome, start_bp, end_bp, output)
    commands <- c(commands, cmd)
  }
  
  # Create script content (removed ibash)
  script_content <- c(
    "#!/bin/bash",
    "module load plink/1.90Beta4.5",
    "",
    commands
  )
  
  # Write script
  writeLines(script_content, "plink_extractions.sh")
  
  # Check if sshpass is available
  sshpass_check <- system("which sshpass", ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (sshpass_check != 0) {
    cat("sshpass not found. Please install it:\n")
    cat("Mac: brew install sshpass\n")
    cat("Linux: sudo apt-get install sshpass\n")
    cat("Or use manual approach: create_plink_script(chunk_files, windows)\n")
    return(NULL)
  }
  
  # Prompt for password
  cat("Enter password for rorytb@ftdcsub.pmacs.upenn.edu: ")
  password <- readline()
  
  # Transfer script
  cat("Transferring script...\n")
  scp_cmd <- sprintf('sshpass -p "%s" scp plink_extractions.sh rorytb@ftdcsub.pmacs.upenn.edu:~/', password)
  scp_result <- system(scp_cmd, ignore.stderr = FALSE)
  
  if (scp_result != 0) {
    cat("Failed to transfer script. Check password and connection.\n")
    return(NULL)
  }
  
  cat("Executing PLINK extractions on LPC...\n")
  # Removed ibash - run directly
  ssh_cmd <- sprintf('sshpass -p "%s" ssh rorytb@ftdcsub.pmacs.upenn.edu "chmod +x plink_extractions.sh && bash plink_extractions.sh"', password)
  ssh_result <- system(ssh_cmd, ignore.stderr = FALSE)
  
  # Cleanup
  unlink("plink_extractions.sh")
  
  if (ssh_result == 0) {
    cat("PLINK extractions completed successfully!\n")
  } else {
    cat("Execution may have failed. Check output above.\n")
    cat("Note: You may need to run this in an interactive session (ibash) manually.\n")
  }
  
  return(merged[, c("CpG_ID", "CpG_Position", "Window_Distance")])
}

# Generate script for manual execution with ibash
create_plink_script <- function(chunk_files, windows_df, script_name = "plink_extractions.sh") {
  
  # Merge to get CpG positions
  merged <- merge(chunk_files, windows_df, by.x = "CpG_ID", by.y = "Input_CpG")
  
  # Generate PLINK commands
  commands <- c()
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
    output <- paste0("/project/ftdc_external/PMBB/genetics/", cpg_id, "_", window_bp, "bp")
    
    cmd <- sprintf('plink --bfile "%s" --chr %s --from-bp %d --to-bp %d --recode A --out "%s"',
                   bfile, chromosome, start_bp, end_bp, output)
    commands <- c(commands, cmd)
  }
  
  # Create script content
  script_content <- c(
    "#!/bin/bash",
    "module load plink/1.90Beta4.5",
    "",
    commands
  )
  
  # Write script
  writeLines(script_content, script_name)
  system(paste("chmod +x", script_name))
  
  cat(paste("Created script:", script_name, "\n"))
  cat("To run on LPC:\n")
  cat(paste("1. scp", script_name, "rorytb@ftdcsub.pmacs.upenn.edu:~/\n"))
  cat("2. ssh rorytb@ftdcsub.pmacs.upenn.edu\n")
  cat("3. ibash\n")
  cat(paste("4. bash", script_name, "\n"))
  
  return(merged[, c("CpG_ID", "CpG_Position", "Window_Distance")])
}

# Usage:
# Automated (tries direct execution): results <- run_plink_extractions(chunk_files, windows)
# Manual (includes ibash step): results <- create_plink_script(chunk_files, windows)