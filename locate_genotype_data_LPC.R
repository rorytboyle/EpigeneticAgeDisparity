# This script locates PMBB genotype data (imputed genotypes as PLINK files) on the LPC
# Author: Rory Boyle rorytboyle@gmail.com
# Date: 17/09/2025

# windows_df is output from get_windows.R

# Function to check if mount point is already mounted
is_mounted <- function(mount_point = "~/lpc_mount") {
  mount_point <- path.expand(mount_point)
  
  # Check using the mount command
  cmd <- "mount"
  mounts <- system(cmd, intern = TRUE)
  
  # Check if our mount point appears in the mount list
  is_mounted <- any(grepl(mount_point, mounts, fixed = TRUE))
  
  if (is_mounted) {
    message("Mount point is currently mounted: ", mount_point)
  } else {
    message("Mount point is not mounted: ", mount_point)
  }
  
  return(is_mounted)
}

# Function to force unmount where mount point is already mounted
force_unmount <- function(mount_point = "~/lpc_mount") {
  mount_point <- path.expand(mount_point)
  
  if (!is_mounted(mount_point)) {
    message("Mount point is not mounted, nothing to unmount")
    return(TRUE)
  }
  
  message("Attempting to force unmount ", mount_point)
  
  os_type <- Sys.info()["sysname"]
  
  if (os_type == "Darwin") {  # macOS
    # Try diskutil first
    cmd <- sprintf("diskutil unmount force %s", mount_point)
    result <- system(cmd, intern = FALSE)
    
    if (result != 0) {
      # Try umount with force flag
      cmd <- sprintf("umount -f %s", mount_point)
      result <- system(cmd, intern = FALSE)
    }
    
    if (result != 0) {
      # Last resort - kill any processes using the mount
      message("Trying to kill processes using the mount point...")
      kill_cmd <- sprintf("lsof +D %s | awk 'NR>1 {print $2}' | xargs kill -9", mount_point)
      system(kill_cmd, intern = FALSE)
      Sys.sleep(2)  # Wait a moment
      
      # Try unmount again
      cmd <- sprintf("diskutil unmount force %s", mount_point)
      result <- system(cmd, intern = FALSE)
    }
  } else {  # Linux
    cmd <- sprintf("umount -f %s", mount_point)
    result <- system(cmd, intern = FALSE)
    
    if (result != 0) {
      cmd <- sprintf("fusermount -uz %s", mount_point)
      result <- system(cmd, intern = FALSE)
    }
  }
  
  if (result == 0) {
    message("Successfully force unmounted ", mount_point)
    return(TRUE)
  } else {
    warning("Failed to force unmount ", mount_point,
            "\nYou may need to run this in terminal: diskutil unmount force ", mount_point)
    return(FALSE)
  }
}

# Function to mount LPC securely using password input 
mount_lpc_secure <- function(username = "rorytb", mount_point = "~/lpc_mount") {
  # Install getPass if not already installed
  if (!requireNamespace("getPass", quietly = TRUE)) {
    install.packages("getPass")
  }
  
  mount_point <- path.expand(mount_point)
  
  # Check if already mounted and try to unmount if necessary
  if (is_mounted(mount_point)) {
    message("Mount point is already in use. Attempting to unmount first...")
    if (!force_unmount(mount_point)) {
      stop("Could not unmount existing mount. Please manually run: diskutil unmount force ", mount_point)
    }
    Sys.sleep(2)  # Wait a moment after unmounting
  }
  
  # Create mount point if it doesn't exist
  if (!dir.exists(mount_point)) {
    dir.create(mount_point, recursive = TRUE)
  }
  
  # Get password securely
  password <- getPass::getPass(paste0("Enter password for ", username, "@ftdcsub.pmacs.upenn.edu: "))
  
  # Use echo to pipe password to sshfs
  cmd <- sprintf("echo '%s' | sshfs -o password_stdin %s@ftdcsub.pmacs.upenn.edu:/static/PMBB/ %s", 
                 password, username, mount_point)
  
  result <- system(cmd, intern = FALSE)
  
  # Clear password from memory
  rm(password)
  
  if (result == 0) {
    message("Successfully mounted LPC at ", mount_point)
    return(TRUE)
  } else {
    warning("Failed to mount LPC")
    return(FALSE)
  }
}

# Function to unmount LPC
unmount_lpc <- function(mount_point = "~/lpc_mount") {
  mount_point <- path.expand(mount_point)
  
  if (!is_mounted(mount_point)) {
    message("Mount point is not mounted, nothing to unmount")
    return(TRUE)
  }
  
  os_type <- Sys.info()["sysname"]
  
  if (os_type == "Darwin") {  # macOS
    # Try diskutil first
    cmd <- sprintf("diskutil unmount %s", mount_point)
    result <- system(cmd, intern = FALSE)
    
    if (result == 0) {
      message("Successfully unmounted ", mount_point)
      return(TRUE)
    } else {
      # Try umount with force flag if diskutil fails
      cmd_alt <- sprintf("umount -f %s", mount_point)
      result_alt <- system(cmd_alt, intern = FALSE)
      
      if (result_alt == 0) {
        message("Successfully unmounted ", mount_point)
        return(TRUE)
      } else {
        warning("Failed to unmount ", mount_point, 
                ". You may need to manually unmount using: diskutil unmount force ", mount_point)
        return(FALSE)
      }
    }
  } else {  # Linux and other systems
    cmd <- sprintf("umount %s", mount_point)
    result <- system(cmd, intern = FALSE)
    
    if (result == 0) {
      message("Successfully unmounted ", mount_point)
      return(TRUE)
    } else {
      cmd_alt <- sprintf("fusermount -u %s", mount_point)
      result_alt <- system(cmd_alt, intern = FALSE)
      
      if (result_alt == 0) {
        message("Successfully unmounted ", mount_point)
        return(TRUE)
      } else {
        warning("Failed to unmount ", mount_point)
        return(FALSE)
      }
    }
  }
}

# Find chunked bed files from windows_df which contains genomic location of CpGs within a given window
find_chunk_files <- function(windows_df, file_path = "~/lpc_mount/PMBB-Release-2024-3.0/Imputed/chunked_bed_files") {
  
  results <- data.frame()
  
  for (i in 1:nrow(windows_df)) {
    cpg_id <- windows_df$Input_CpG[i]
    chromosome <- windows_df$Chromosome[i]
    start_pos <- windows_df$Window_Start[i]
    end_pos <- windows_df$Window_End[i]
    
    # Find ONLY .bed files for this chromosome (fixed regex pattern)
    files <- list.files(file_path, pattern = paste0(".*", chromosome, ".*\\.bed$"), full.names = TRUE)
    
    # Check each .bed file for overlap
    for (file in files) {
      filename <- basename(file)
      
      # Extract chunk coordinates: chr*_chunk*_start_end.bed
      if (grepl(paste0(chromosome, "_chunk[0-9]+_([0-9]+)_([0-9]+)\\.bed$"), filename)) {
        matches <- regmatches(filename, regexec(paste0(chromosome, "_chunk[0-9]+_([0-9]+)_([0-9]+)\\.bed$"), filename))
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

# wrapper function
process_chunk_files <- function(windows_df, 
                                username = "rorytb",
                                mount_point = "~/lpc_mount",
                                file_path = "~/lpc_mount/PMBB-Release-2024-3.0/Imputed/chunked_bed_files",
                                use_secure = TRUE) {
  
  # Mount LPC (will handle existing mounts automatically)
  if (use_secure) {
    mount_success <- mount_lpc_secure(username, mount_point)
  } else {
    stop("Please use use_secure = TRUE for password-based mounting")
  }
  
  if (!mount_success) {
    stop("Failed to mount LPC. Exiting.")
  }
  
  # Use tryCatch to ensure unmounting even if there's an error
  chunk_files <- tryCatch({
    message("Processing chunk files...")
    find_chunk_files(windows_df, file_path)
  }, error = function(e) {
    message("Error occurred: ", e$message)
    NULL
  }, finally = {
    # Always try to unmount
    unmount_lpc(mount_point)
  })
  
  return(chunk_files)
}

# Usage ####
# Clean up any existing mount
force_unmount("~/lpc_mount")

# find chunked files
chunk_files <- process_chunk_files(windows, use_secure = TRUE)