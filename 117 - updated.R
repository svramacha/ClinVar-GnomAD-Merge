library(httr)
library(jsonlite)
library(VariantAnnotation)

# Edited query code

# Checked what can be queried via terminal
query <- '{"query": "{ datasets { name description } }"}'

# Check if im getting response from query
url <- "https://gnomad.broadinstitute.org/api" 
response <- POST(url, body = query, encode = "raw", content_type_json())
content <- content(response, "parsed")
print(content)

# Define regions 
regions <- data.frame(
  chrom = "10",
  start = 87863114,
  stop = 87971941
)

print(regions)

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





# Run this one!!
if(!require(conflicted)){install.packages("conflicted")}
library(conflicted)
library(httr)
library(jsonlite)
library(VariantAnnotation)

conflicts_prefer(httr::content) # Set preference for httr::content

# Edited query code

# Check what can be queried via terminal
query <- '{"query": "{ datasets { name description } }"}'

# Check if I'm getting response from query
url <- "https://gnomad.broadinstitute.org/api"
response <- POST(url, body = query, encode = "raw", content_type_json())

# **Correct usage of content() to get the JSON text:**
content <- content(response, "text")  # Get the raw JSON text

# **Parse the JSON text using jsonlite (assuming it's valid JSON):**
data <- fromJSON(content)  # Parse the JSON into an R object (list)

# Print the data (for initial testing)
print(data)

# Define regions
regions <- data.frame(
  chrom = "10",
  start = 87863114,
  stop = 87971941
)

print(regions)

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

# Run the API

response <- POST(url, 
                 body = body, 
                 encode = "raw", 
                 content_type("application/json"), 
                 timeout(60))

if (status_code(response) == 200) {
  content_text <- httr::content(response, "text")
  data <- fromJSON(content_text)
  
  variants_df <- data$data$region$variants
  
  if (nrow(variants_df) > 0) {
    # Correctly handle list columns and create new columns
    variants_df$flags <- sapply(variants_df$flags, function(x) paste(x, collapse = ","))
    variants_df$filters <- sapply(variants_df$joint$filters, function(x) paste(x, collapse = ","))
    variants_df$ac <- variants_df$joint$ac
    variants_df$an <- variants_df$joint$an
    variants_df$faf95_max <- variants_df$joint$fafmax$faf95_max
    variants_df$faf95_max_gen_anc <- variants_df$joint$fafmax$faf95_max_gen_anc
    variants_df$joint <- NULL
    
    # Add region information
    variants_df$chrom <- regions$chrom[1]
    variants_df$start <- regions$start[1]
    variants_df$stop <- regions$stop[1]
    
    print("Data frame created successfully:")
    print(head(variants_df))
    print(dim(variants_df)) # Print dimensions before writing to CSV
    
    write.csv(variants_df, "variants_data.csv", row.names = FALSE)
  } else {
    print("No variants found for this region (empty data frame).")
  }
} else {
  print(paste("Error:", status_code(response)))
  error_content <- httr::content(response, as = "text")
  print(error_content)
}


