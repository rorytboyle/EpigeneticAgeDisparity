# This script locates PMBB genotype data (imputed genotypes as PLINK files) on the LPC
# Author: Rory Boyle rorytboyle@gmail.com
# Date: 09/09/2025

# windows_df is output from get_windows.R

# Try to mount PMBB

# # Simple function to find overlapping chunk files
find_chunk_files <- function(windows_df, file_path = "~/lpc_mount/PMBB-Release-2024-3.0/Imputed/chunked_bed_files/") {
  
  results <- data.frame()
  
  for (i in 1:nrow(windows_df)) {
    cpg_id <- windows_df$Input_CpG[i]
    chromosome <- windows_df$Chromosome[i]
    start_pos <- windows_df$Window_Start[i]
    end_pos <- windows_df$Window_End[i]
    
    # Find files for this chromosome
    files <- list.files(file_path, pattern = paste0("*", chromosome, "*"), full.names = TRUE)
    
    # Check each file for overlap
    for (file in files) {
      filename <- basename(file)
      
      # Extract chunk coordinates: chr*_chunk*_start_end
      if (grepl(paste0(chromosome, "_chunk[0-9]+_([0-9]+)_([0-9]+)"), filename)) {
        matches <- regmatches(filename, regexec(paste0(chromosome, "_chunk[0-9]+_([0-9]+)_([0-9]+)"), filename))
        chunk_start <- as.numeric(matches[[1]][2])
        chunk_end <- as.numeric(matches[[1]][3])
        
        # Check overlap
        if (chunk_start <= end_pos && chunk_end >= start_pos) {
          results <- rbind(results, data.frame(
            CpG_ID = cpg_id,
            File_Path = file,
            Chunk_Start = chunk_start,
            Chunk_End = chunk_end
          ))
        }
      }
    }
  }
  
  return(results)
}

# Example usage ####

# Mount LPC
# sshfs rorytb@ftdcsub.pmacs.upenn.edu:/static/PMBB/ ~/lpc_mount       

# Run function
chunk_files <- find_chunk_files(windows)

# Unmount LPC
# umount ~/lpc_mount
