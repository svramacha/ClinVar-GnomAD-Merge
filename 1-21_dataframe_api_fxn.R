# Load necessary libraries
library(GenomicRanges)
library(httr)
library(jsonlite)
library(Rsamtools)
library(VariantAnnotation)
if(!require(conflicted)){install.packages("conflicted")}
library(conflicted)

conflicts_prefer(httr::content)





library(httr)
library(jsonlite)

# Function to query the gnomAD API, ClinVar API, and CADD API
fetch_variant_data <- function(chrom, start, stop) {
  
  # Define gnomAD query
  gnomad_query <- '
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
  
  # Define query variables for gnomAD
  gnomad_variables <- list(
    chrom = chrom,
    start = start,
    stop = stop
  )
  
  # Create the POST request body for gnomAD
  gnomad_body <- toJSON(list(
    query = gnomad_query,
    variables = gnomad_variables
  ), auto_unbox = TRUE)
  
  # Define gnomAD API URL
  gnomad_url <- "https://gnomad.broadinstitute.org/api"
  
  # Send POST request to gnomAD
  response_gnomad <- POST(gnomad_url, 
                          body = gnomad_body, 
                          encode = "raw", 
                          content_type_json(), 
                          timeout(90))
  
  # Check if the response is successful
  if (status_code(response_gnomad) == 200) {
    content_gnomad <- content(response_gnomad, "text")
    data_gnomad <- fromJSON(content_gnomad)
    
    variants_df <- data_gnomad$data$region$variants
    
    # Handle missing variants
    if (length(variants_df) == 0) {
      print("No variants found for this region.")
      return(NULL)
    }
    
    # Clean up gnomAD data (processing flags, filters, and adding region info)
    variants_df$flags <- sapply(variants_df$flags, function(x) paste(x, collapse = ","))
    variants_df$filters <- sapply(variants_df$joint$filters, function(x) paste(x, collapse = ","))
    variants_df$ac <- variants_df$joint$ac
    variants_df$an <- variants_df$joint$an
    variants_df$faf95_max <- variants_df$joint$fafmax$faf95_max
    variants_df$faf95_max_gen_anc <- variants_df$joint$fafmax$faf95_max_gen_anc
    variants_df$joint <- NULL
    
    variants_df$chrom <- chrom
    variants_df$start <- start
    variants_df$stop <- stop
    
    # Now fetch ClinVar and CADD annotations for each variant in gnomAD
    clinvar_data <- sapply(variants_df$variant_id, function(variant_id) {
      clinvar_url <- paste0("https://api.clinvar.gov/v1/variant/", variant_id)
      clinvar_response <- GET(clinvar_url)
      
      if (status_code(clinvar_response) == 200) {
        clinvar_content <- content(clinvar_response, "text")
        clinvar_json <- fromJSON(clinvar_content)
        return(clinvar_json)
      } else {
        return(NULL)
      }
    })
    
    cadd_data <- sapply(variants_df$variant_id, function(variant_id) {
      cadd_url <- paste0("https://cadd.gs.washington.edu/api/v1/variant/", variant_id)
      cadd_response <- GET(cadd_url)
      
      if (status_code(cadd_response) == 200) {
        cadd_content <- content(cadd_response, "text")
        cadd_json <- fromJSON(cadd_content)
        return(cadd_json)
      } else {
        return(NULL)
      }
    })
    
    # Merge ClinVar and CADD data into the variants data frame
    variants_df$clinvar_annotation <- clinvar_data
    variants_df$cadd_score <- sapply(cadd_data, function(x) x$score)
    
    # Print results
    print("Data frame created successfully:")
    print(head(variants_df))
    
    # Save to CSV
    write.csv(variants_df, "variants_with_annotations.csv", row.names = FALSE)
    
  } else {
    print(paste("Error fetching data from gnomAD. Status code:", status_code(response_gnomad)))
    error_content <- content(response_gnomad, "text")
    print(error_content)
  }
}

# Example usage of the function
fetch_variant_data(chrom = "10", start = 87863114, stop = 87971941)
