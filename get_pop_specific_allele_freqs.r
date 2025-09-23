# This script executes gnomAD API queries to get population-specific allele frequencies for variants
# around each cpg site in windows_df.
# Author: Rory Boyle & Nadia Dehghani rorytboyle@gmail.com
# Date: 23/09/2025

# windows_df is output from get_windows.R

library(httr2)
library(jsonlite)
library(tidyverse)

# Function to query gnomAD API for multiple genomic regions
query_gnomad_regions <- function(windows_df, 
                                 dataset = "gnomad_r4", 
                                 reference_genome = "GRCh38",
                                 save_individual_files = TRUE,
                                 output_dir = "gnomad_queries",
                                 delay_seconds = 0.5) {
  
  # Create output directory if it doesn't exist
  if (save_individual_files && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Initialize results list
  all_results <- list()
  failed_queries <- c()
  
  cat(sprintf("Querying gnomAD for %d regions...\n", nrow(windows_df)))
  
  # Process each row in windows_df
  for (i in 1:nrow(windows_df)) {
    
    # Extract information for this region
    cpg_id <- windows_df$Input_CpG[i]
    chrom <- gsub("chr", "", windows_df$Chromosome[i])  # Remove chr prefix
    start_pos <- windows_df$Window_Start[i]
    end_pos <- windows_df$Window_End[i]
    
    cat(sprintf("Processing %d/%d: %s (chr%s:%d-%d)...\n", 
                i, nrow(windows_df), cpg_id, chrom, start_pos, end_pos))
    
    # Construct GraphQL query
    graphql_query <- sprintf('
    query getVariantsInRegion {
      region(reference_genome: %s, chrom: "%s", start: %d, stop: %d) {
        variants (dataset: %s) {
          variant_id
          pos
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
    
    # Execute query
    tryCatch({
      
      # Make API request
      response <- request("https://gnomad.broadinstitute.org/api") %>%
        req_headers("Content-Type" = "application/json") %>%
        req_body_json(list(query = graphql_query)) %>%
        req_perform()
      
      # Parse response
      response_data <- response %>% resp_body_json()
      
      # Check for errors
      if (!is.null(response_data$errors)) {
        warning(sprintf("GraphQL errors for %s: %s", cpg_id, 
                        paste(sapply(response_data$errors, function(x) x$message), collapse = "; ")))
        failed_queries <- c(failed_queries, cpg_id)
        next
      }
      
      # Process the data
      if (!is.null(response_data$data$region$variants) && 
          length(response_data$data$region$variants) > 0) {
        
        # Convert to tidy format
        tidy_data <- process_gnomad_response(response_data, cpg_id, chrom, start_pos, end_pos)
        
        # Save individual file if requested
        if (save_individual_files) {
          filename <- file.path(output_dir, paste0(cpg_id, "_gnomad.csv"))
          write_csv(tidy_data, filename)
          cat(sprintf("  Saved %d variants to %s\n", nrow(tidy_data), filename))
        }
        
        # Add to results list
        all_results[[cpg_id]] <- tidy_data
        
      } else {
        cat(sprintf("  No variants found for %s\n", cpg_id))
        # Still add empty dataframe to maintain structure
        empty_df <- create_empty_gnomad_df(cpg_id, chrom, start_pos, end_pos)
        all_results[[cpg_id]] <- empty_df
        
        if (save_individual_files) {
          filename <- file.path(output_dir, paste0(cpg_id, "_gnomad.csv"))
          write_csv(empty_df, filename)
        }
      }
      
      # Add delay to be respectful to the API
      if (delay_seconds > 0 && i < nrow(windows_df)) {
        Sys.sleep(delay_seconds)
      }
      
    }, error = function(e) {
      warning(sprintf("Error querying %s: %s", cpg_id, e$message))
      failed_queries <- c(failed_queries, cpg_id)
    })
  }
  
  # Summary
  cat(sprintf("\nCompleted queries for %d regions\n", nrow(windows_df)))
  cat(sprintf("Successful: %d\n", length(all_results)))
  cat(sprintf("Failed: %d\n", length(failed_queries)))
  
  if (length(failed_queries) > 0) {
    cat("Failed queries:", paste(failed_queries, collapse = ", "), "\n")
  }
  
  # Return results
  return(list(
    data = all_results,
    failed_queries = failed_queries,
    summary = data.frame(
      total_regions = nrow(windows_df),
      successful = length(all_results),
      failed = length(failed_queries),
      output_directory = ifelse(save_individual_files, output_dir, NA)
    )
  ))
}

# Function to process gnomAD API response into tidy format
process_gnomad_response <- function(response_data, cpg_id, chrom, start_pos, end_pos) {
  
  variants <- response_data$data$region$variants
  
  if (length(variants) == 0) {
    return(create_empty_gnomad_df(cpg_id, chrom, start_pos, end_pos))
  }
  
  # Convert to dataframe
  variant_df <- map_dfr(variants, function(variant) {
    
    # Extract basic variant info
    basic_info <- data.frame(
      region_id = cpg_id,
      chromosome = chrom,
      region_start = start_pos,
      region_end = end_pos,
      variant_id = variant$variant_id %||% NA,
      position = variant$pos %||% NA,
      rsids = ifelse(is.null(variant$rsids) || length(variant$rsids) == 0, 
                     NA, paste(variant$rsids, collapse = ";")),
      stringsAsFactors = FALSE
    )
    
    # Extract population data
    if (!is.null(variant$genome) && !is.null(variant$genome$populations)) {
      pop_data <- map_dfr(variant$genome$populations, function(pop) {
        data.frame(
          population_id = pop$id %||% NA,
          ac = pop$ac %||% NA,
          an = pop$an %||% NA,
          stringsAsFactors = FALSE
        )
      })
      
      # Reshape population data to wide format
      pop_wide <- pop_data %>%
        pivot_wider(
          names_from = population_id,
          values_from = c(ac, an),
          names_sep = "_"
        )
      
      # Combine basic info with population data
      result <- bind_cols(basic_info, pop_wide)
      
    } else {
      result <- basic_info
    }
    
    return(result)
  })
  
  # Calculate allele frequencies
  pop_columns <- names(variant_df)[grepl("^ac_", names(variant_df))]
  
  for (pop_col in pop_columns) {
    pop_name <- gsub("^ac_", "", pop_col)
    an_col <- paste0("an_", pop_name)
    af_col <- paste0("af_", pop_name)
    
    if (an_col %in% names(variant_df)) {
      variant_df[[af_col]] <- ifelse(variant_df[[an_col]] > 0, 
                                     variant_df[[pop_col]] / variant_df[[an_col]], 
                                     NA)
    }
  }
  
  return(variant_df)
}

# Helper function to create empty dataframe with consistent structure
create_empty_gnomad_df <- function(cpg_id, chrom, start_pos, end_pos) {
  data.frame(
    region_id = cpg_id,
    chromosome = chrom,
    region_start = start_pos,
    region_end = end_pos,
    variant_id = character(0),
    position = numeric(0),
    rsids = character(0),
    stringsAsFactors = FALSE
  )
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Function to combine all results into a single dataframe
combine_gnomad_results <- function(gnomad_results) {
  if (length(gnomad_results$data) == 0) {
    warning("No successful queries to combine")
    return(NULL)
  }
  
  # Combine all dataframes
  combined_df <- map_dfr(gnomad_results$data, ~ .x, .id = "query_name")
  
  return(combined_df)
}

# Function to get summary statistics
summarize_gnomad_results <- function(gnomad_results) {
  if (length(gnomad_results$data) == 0) {
    return(data.frame(region_id = character(0), variant_count = numeric(0)))
  }
  
  summary_stats <- map_dfr(gnomad_results$data, function(df) {
    data.frame(
      region_id = df$region_id[1],
      variant_count = nrow(df),
      stringsAsFactors = FALSE
    )
  })
  
  return(summary_stats)
}

# Usage: ####

# Run queries for all regions in windows_df
gnomad_results <- query_gnomad_regions(windows_df,
                                      save_individual_files = TRUE,
                                      output_dir = "/Users/rorytb/Library/CloudStorage/Box-Box/PennMedicineBiobank/DNAmethylation/gnomad_output",
                                      delay_seconds = 0.5)
 
# View summary
print(gnomad_results$summary)

# Get summary statistics
variant_summary <- summarize_gnomad_results(gnomad_results)
print(variant_summary)

# Combine all results into single dataframe
combined_data <- combine_gnomad_results(gnomad_results)
 
# Access individual cpg data
cg25243766_data <- gnomad_results$data$cg25243766