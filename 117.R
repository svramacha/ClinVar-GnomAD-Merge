# Load necessary libraries
if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("VariantAnnotation", quietly = TRUE)) BiocManager::install("VariantAnnotation")
if (!requireNamespace("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")
if (!requireNamespace("Rsamtools", quietly = TRUE)) BiocManager::install("Rsamtools")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(httr)
library(jsonlite)
library(VariantAnnotation)
library(GenomicRanges)
library(Rsamtools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)

# Step 1: Set up the API endpoint
url <- "https://gnomad.broadinstitute.org/api"

# Step 2: Define the GraphQL query
query <- '
query GeneVariantsFilteringAF($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      ref
      alt
      flags
      joint {
        ac
        an
        filters
        fafmax {
          faf95_max
          faf95_max_gen_anc
        }
      }
    }
  }
}
'

# Step 3: Define regions for CHEK2
regions <- data.frame(
  chrom = "22",
  start = 28687820,
  stop = 28742014
)

print(regions)

# Step 4: Initialize storage for results
all_results <- list()

# Step 5: Iterate over each region with retry logic
max_retries <- 3
for (i in seq_len(nrow(regions))) {
  variables <- list(
    chrom = regions$chrom[i],
    start = regions$start[i],
    stop = regions$stop[i]
  )
  
  # Create the POST request body
  body <- toJSON(list(
    query = query,
    variables = variables
  ), auto_unbox = TRUE)
  
  for (attempt in 1:max_retries) {
    response <- tryCatch({
      POST(
        url,
        body = body,
        encode = "raw",  # Use raw since the body is already JSON
        content_type("application/json"),
        timeout(60)  # Set timeout to 60 seconds
      )
    }, error = function(e) NULL)
    
    if (!is.null(response) && status_code(response) == 200) {
      data <- content(response, "parsed")
      variants <- data$data$region$variants
      if (!is.null(variants)) {
        # Convert to data frame and store
        results_df <- do.call(rbind, lapply(variants, as.data.frame))
        results_df$chrom <- regions$chrom[i]
        results_df$start <- regions$start[i]
        results_df$stop <- regions$stop[i]
        all_results[[i]] <- results_df
      } else {
        print(paste("No variants found for region", i))
      }
      break
    } else {
      if (attempt == max_retries) {
        print(paste("Error with API call for region", i, ": HTTP status", status_code(response)))
      } else {
        Sys.sleep(2)  # Wait before retrying
      }
    }
  }
}

# Step 6: Combine results into a single data frame
if (length(all_results) > 0) {
  combined_results <- do.call(rbind, all_results)
  print(head(combined_results))
} else {
  print("No data retrieved for any region.")
}

# Step 7: Save results to a CSV file
if (length(all_results) > 0) {
  write.csv(combined_results, "gnomad_results_HBOC.csv", row.names = FALSE)
}




# Edited query code

# Checked what can be queried via terminal
query <- '{"query": "{ datasets { name description } }"}'

# Check if im getting response from query
response <- POST(url, body = query, encode = "raw", content_type_json())
content <- content(response, "parsed")
print(content)
print(response)

# Query with API and save into dataframe 
query <- '
query GeneVariantsFilteringAF($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      ref
      alt
      flags
      joint {
        ac
        an
        filters
        fafmax {
          faf95_max
          faf95_max_gen_anc
        }
      }
    }
  }
}
'
# Query variables
variables <- list(
  chrom = regions$chrom[1],
  start = regions$start[1],
  stop = regions$stop[1]
)

# Create the POST request body
body <- toJSON(list(
  query = query,
  variables = variables
), auto_unbox = TRUE)

# Make the API call
response <- POST(
  url,
  body = body,
  encode = "raw", 
  content_type("application/json"),
  timeout(60)
)

# Check the response
if (status_code(response) == 200) {
  data <- content(response, "parsed")
  variants <- data$data$region$variants
  
  if (!is.null(variants)) {
    # Initialize a list to store data frames for each variant
    variant_list <- list()
    
    # Iterate over each variant to create a consistent data frame
    for (variant in variants) {
      # Safely extract each field, using NA if the field is missing
      variant_data <- data.frame(
        variant_id = ifelse(!is.null(variant$variant_id), variant$variant_id, NA),
        ref = ifelse(!is.null(variant$ref), variant$ref, NA),
        alt = ifelse(!is.null(variant$alt), variant$alt, NA),
        flags = ifelse(!is.null(variant$flags), paste(variant$flags, collapse = ","), NA),
        ac = ifelse(!is.null(variant$joint$ac), variant$joint$ac, NA),
        an = ifelse(!is.null(variant$joint$an), variant$joint$an, NA),
        filters = ifelse(!is.null(variant$joint$filters), paste(variant$joint$filters, collapse = ","), NA),
        faf95_max = ifelse(!is.null(variant$joint$fafmax$faf95_max), variant$joint$fafmax$faf95_max, NA),
        faf95_max_gen_anc = ifelse(!is.null(variant$joint$fafmax$faf95_max_gen_anc), variant$joint$fafmax$faf95_max_gen_anc, NA)
      )
      
      # Add the variant data to the list
      variant_list <- append(variant_list, list(variant_data))
    }
    
    # Combine all variant data frames into a single data frame
    variants_df <- do.call(rbind, variant_list)
    
    # Add region information to the data frame
    variants_df$chrom <- regions$chrom[1]
    variants_df$start <- regions$start[1]
    variants_df$stop <- regions$stop[1]
    
    # Print the first few rows of the data frame
    print("Data frame created successfully:")
    print(head(variants_df))
    
    # Optionally, save the data frame to a CSV file
    write.csv(variants_df, "variants_data.csv", row.names = FALSE)
  } else {
    print("No variants found for this region.")
  }
} else {
  print(paste("Error: ", status_code(response)))
}
