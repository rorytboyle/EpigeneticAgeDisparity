# This script executes gnomAD API queries to get population-specific allele frequencies for variants
# around each cpg site in windows_df.
# Author: Rory Boyle & Nadia Dehghani rorytboyle@gmail.com
# Date: 23/09/2025

# windows_df is output from get_windows.R

# Create GraphQL query for gnomAD API
create_gnomad_query <- function(chromosome, start_pos, end_pos, 
                                dataset = "gnomad_r4", 
                                reference_genome = "GRCh38") {
  
  chrom <- gsub("chr", "", chromosome)
  
  query <- sprintf('
  query getVariantsInRegion {
    region(reference_genome: %s, chrom: "%s", start: %d, stop: %d) {
      variants (dataset: %s) {
        variant_id
        pos
        ref
        alt
        rsids
        genome {
          populations {
            id
            ac
            an
          }
        }
      }
    }
  }', reference_genome, chrom, start_pos, end_pos, dataset)
  
  return(query)
}

# Query gnomAD for all population frequencies
query_gnomad_all_populations <- function(windows_df, 
                                         dataset = "gnomad_r4",
                                         save_individual_files = TRUE,
                                         output_dir = "snp_lists",
                                         delay_seconds = 0.5) {
  
  library(httr2)
  library(jsonlite)
  library(tidyverse)
  
  if (save_individual_files && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  all_snp_lists <- list()
  
  cat(sprintf("Querying gnomAD for %d regions...\n", nrow(windows_df)))
  
  for (i in 1:nrow(windows_df)) {
    cpg_id <- windows_df$Input_CpG[i]
    chrom <- gsub("chr", "", windows_df$Chromosome[i])
    start_pos <- windows_df$Window_Start[i]
    end_pos <- windows_df$Window_End[i]
    
    cat(sprintf("Processing %d/%d: %s (chr%s:%d-%d)...\n", 
                i, nrow(windows_df), cpg_id, chrom, start_pos, end_pos))
    
    tryCatch({
      
      query <- create_gnomad_query(chrom, start_pos, end_pos, dataset)
      
      response <- request("https://gnomad.broadinstitute.org/api") %>%
        req_headers("Content-Type" = "application/json") %>%
        req_body_json(list(query = query)) %>%
        req_perform()
      
      response_data <- response %>% resp_body_json()
      
      if (!is.null(response_data$errors)) {
        warning(sprintf("GraphQL errors for %s: %s", cpg_id, 
                        paste(sapply(response_data$errors, function(x) x$message), collapse = "; ")))
        next
      }
      
      if (!is.null(response_data$data$region$variants) && 
          length(response_data$data$region$variants) > 0) {
        
        snp_list <- process_gnomad_response(response_data, cpg_id)
        
        if (save_individual_files && !is.null(snp_list)) {
          filename <- file.path(output_dir, paste0(cpg_id, "_gnomAD_af.csv"))
          write_csv(snp_list, filename)
          cat(sprintf("  Saved %d SNPs with %d populations to %s\n", 
                      nrow(snp_list), 
                      sum(grepl("^af_", names(snp_list))), 
                      filename))
        }
        
        all_snp_lists[[cpg_id]] <- snp_list
        
      } else {
        cat(sprintf("  No variants found for %s\n", cpg_id))
        all_snp_lists[[cpg_id]] <- create_empty_snp_df(cpg_id)
      }
      
      if (delay_seconds > 0 && i < nrow(windows_df)) {
        Sys.sleep(delay_seconds)
      }
      
    }, error = function(e) {
      warning(sprintf("Error processing %s: %s", cpg_id, e$message))
    })
  }
  
  cat(sprintf("\nCompleted queries for %d regions\n", nrow(windows_df)))
  cat(sprintf("Successful: %d\n", length(all_snp_lists)))
  
  return(all_snp_lists)
}

# Process gnomAD response to include all population frequencies
process_gnomad_response <- function(response_data, cpg_id) {
  
  variants <- response_data$data$region$variants
  
  if (length(variants) == 0) {
    return(create_empty_snp_df(cpg_id))
  }
  
  # Get all unique population IDs
  all_populations <- unique(unlist(lapply(variants, function(v) {
    if (!is.null(v$genome$populations)) {
      sapply(v$genome$populations, function(p) p$id %||% "unknown")
    }
  })))
  
  # Sort alphabetically but keep "remaining" at the end
  remaining_pops <- all_populations[all_populations == "remaining"]
  other_pops <- sort(all_populations[all_populations != "remaining"])
  all_populations <- c(other_pops, remaining_pops)
  
  cat(sprintf("    Found populations: %s\n", paste(all_populations, collapse = ", ")))
  
  # Convert to SNP dataframe
  snp_df <- map_dfr(variants, function(variant) {
    
    # Extract chromosome from variant_id
    variant_parts <- strsplit(variant$variant_id, "-")[[1]]
    chromosome <- ifelse(length(variant_parts) >= 4, variant_parts[1], NA)
    
    # Basic SNP info
    snp_row <- data.frame(
      region_id = cpg_id,
      chr = chromosome,
      position = variant$pos %||% NA,
      allele1 = variant$ref %||% NA,
      allele2 = variant$alt %||% NA,
      rsid = ifelse(is.null(variant$rsids) || length(variant$rsids) == 0, 
                    NA, variant$rsids[[1]]),
      stringsAsFactors = FALSE
    )
    
    # Initialize population columns
    for (pop in all_populations) {
      snp_row[[paste0("ac_", pop)]] <- NA_integer_
      snp_row[[paste0("an_", pop)]] <- NA_integer_
      snp_row[[paste0("af_", pop)]] <- NA_real_
    }
    
    # Add population data
    if (!is.null(variant$genome$populations)) {
      for (pop_data in variant$genome$populations) {
        pop_id <- pop_data$id %||% "unknown"
        ac <- pop_data$ac %||% 0
        an <- pop_data$an %||% 0
        af <- ifelse(an > 0, ac / an, NA)
        
        snp_row[[paste0("ac_", pop_id)]] <- ac
        snp_row[[paste0("an_", pop_id)]] <- an
        snp_row[[paste0("af_", pop_id)]] <- af
      }
    }
    
    return(snp_row)
  })
  
  # Clean and sort
  snp_df <- snp_df %>%
    filter(!is.na(position)) %>%
    arrange(as.numeric(chr), position)
  
  return(snp_df)
}

# Create empty SNP dataframe
create_empty_snp_df <- function(cpg_id) {
  data.frame(
    region_id = cpg_id,
    chr = character(0),
    position = numeric(0),
    allele1 = character(0),
    allele2 = character(0),
    rsid = character(0),
    stringsAsFactors = FALSE
  )
}

# Combine all SNP lists into master list
combine_snp_lists <- function(snp_lists) {
  
  if (length(snp_lists) == 0) {
    warning("No SNP lists to combine")
    return(NULL)
  }
  
  cat("Combining SNP lists...\n")
  
  # Get all column names
  all_columns <- unique(unlist(lapply(snp_lists, names)))
  
  # Standardize all dataframes to have same columns
  snp_lists_std <- lapply(snp_lists, function(df) {
    missing_cols <- setdiff(all_columns, names(df))
    
    if (length(missing_cols) > 0) {
      for (col in missing_cols) {
        if (grepl("^(af_|ac_|an_)", col)) {
          df[[col]] <- NA_real_
        } else if (col == "position") {
          df[[col]] <- NA_integer_
        } else {
          df[[col]] <- NA_character_
        }
      }
    }
    
    return(df[, all_columns])
  })
  
  # Combine dataframes
  combined <- map_dfr(snp_lists_std, ~ .x, .id = "source_region")
  
  if (nrow(combined) > 0) {
    # Remove duplicates
    combined <- combined %>%
      distinct(chr, position, allele1, allele2, .keep_all = TRUE) %>%
      arrange(as.numeric(chr), position)
    
    cat(sprintf("Combined %d unique SNPs\n", nrow(combined)))
    
    # Report populations (in sorted order)
    all_af_cols <- names(combined)[grepl("^af_", names(combined))]
    populations <- gsub("^af_", "", all_af_cols)
    
    # Sort alphabetically with "remaining" at end
    remaining_pops <- populations[populations == "remaining"]
    other_pops <- sort(populations[populations != "remaining"])
    sorted_populations <- c(other_pops, remaining_pops)
    
    cat(sprintf("Populations: %s\n", paste(sorted_populations, collapse = ", ")))
  }
  
  return(combined)
}

# Summarize population data completeness
summarize_populations <- function(combined_snps) {
  
  if (is.null(combined_snps) || nrow(combined_snps) == 0) {
    cat("No data to summarize\n")
    return(NULL)
  }
  
  # Get populations (sorted alphabetically with "remaining" at end)
  all_populations <- gsub("^af_", "", names(combined_snps)[grepl("^af_", names(combined_snps))])
  remaining_pops <- all_populations[all_populations == "remaining"]
  other_pops <- sort(all_populations[all_populations != "remaining"])
  populations <- c(other_pops, remaining_pops)
  
  # Summary for each population
  summary_df <- map_dfr(populations, function(pop) {
    af_col <- paste0("af_", pop)
    an_col <- paste0("an_", pop)
    
    af_values <- combined_snps[[af_col]]
    an_values <- if (an_col %in% names(combined_snps)) combined_snps[[an_col]] else rep(NA, length(af_values))
    
    data.frame(
      population = pop,
      total_snps = nrow(combined_snps),
      snps_with_data = sum(!is.na(af_values)),
      snps_with_freq = sum(af_values > 0, na.rm = TRUE),
      mean_af = round(mean(af_values, na.rm = TRUE), 4),
      median_af = round(median(af_values, na.rm = TRUE), 4),
      mean_sample_size = round(mean(an_values, na.rm = TRUE), 0),
      completeness_pct = round(100 * sum(!is.na(af_values)) / nrow(combined_snps), 1),
      stringsAsFactors = FALSE
    )
  })
  
  print(summary_df)
  
  return(summary_df)
}

# Helper function
`%||%` <- function(x, y) if (is.null(x)) y else x

# Usage:

# Query gnomAD for all populations
results <- query_gnomad_all_populations(windows, output_dir = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/gnomad_output")

# Combine into master list
master_snps <- combine_snp_lists(results)

# Get population summary
pop_summary <- summarize_populations(master_snps)

# View results
head(master_snps)
