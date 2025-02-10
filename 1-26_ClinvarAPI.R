library(httr)
library(jsonlite)
library(dplyr)
library(data.table)
library(XML)

library(httr)
library(jsonlite)
library(XML)
library(dplyr)



get_clinvar_variants_by_gene <- function(gene_name, desired_columns = NULL) {
  library(httr)
  library(jsonlite)
  library(XML)
  
  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  
  # First, get the total count of variants for the gene
  esearch_count_url <- paste0(base_url, "esearch.fcgi?db=clinvar&term=", gene_name, "[gene]&rettype=count&retmode=json")
  count_response <- GET(esearch_count_url)
  count_json <- fromJSON(content(count_response, "text"))
  total_records <- as.numeric(count_json$esearchresult$count)
  
  if (total_records == 0) {
    print(paste("No ClinVar records found for gene:", gene_name))
    return(NULL)
  }
  
  print(paste("Total ClinVar records for", gene_name, ":", total_records))
  
  all_ids <- c()
  retmax <- 500  # Max records per request
  
  # Loop through pages of results
  for (retstart in seq(0, total_records, by = retmax)) {
    esearch_url <- paste0(base_url, "esearch.fcgi?db=clinvar&term=", gene_name, "[gene]&retmax=", retmax, "&retstart=", retstart)
    esearch_response <- GET(esearch_url)
    
    if (status_code(esearch_response) != 200) {
      print(paste("Error in esearch request:", status_code(esearch_response)))
      next
    }
    
    esearch_content <- content(esearch_response, "text")
    esearch_xml <- xmlParse(esearch_content)
    ids <- xpathSApply(esearch_xml, "//IdList/Id", xmlValue)
    
    all_ids <- c(all_ids, ids)
    
    Sys.sleep(1)  # Add a 1-second delay after each request
  }
  
  if (length(all_ids) == 0) {
    print(paste("No ClinVar records found for gene:", gene_name))
    return(NULL)
  }
  
  # Retrieve summary details for all IDs in batches
  batch_size <- 100
  filtered_results <- list()
  
  for (batch_start in seq(1, length(all_ids), by = batch_size)) {
    batch_ids <- all_ids[batch_start:min(batch_start + batch_size - 1, length(all_ids))]
    esummary_url <- paste0(base_url, "esummary.fcgi?db=clinvar&id=", paste(batch_ids, collapse = ","), "&retmode=json")
    
    esummary_response <- GET(esummary_url)
    if (status_code(esummary_response) != 200) {
      print(paste("Error in esummary request:", status_code(esummary_response)))
      next
    }
    
    esummary_json <- fromJSON(content(esummary_response, "text"))
    results_list <- esummary_json$result[-1]
    
    for (x in results_list) {
      if ("obj_type" %in% names(x) && "molecular_consequence_list" %in% names(x)) {
        obj_type <- x[["obj_type"]]
        molecular_consequence_list <- x[["molecular_consequence_list"]]
        
        if (is.list(molecular_consequence_list)) {
          molecular_consequence_list <- unlist(molecular_consequence_list)
        }
        
        if (obj_type == "single nucleotide variant" && any(molecular_consequence_list == "missense variant")) {
          tryCatch({
            extracted_data <- list()
            if (is.null(desired_columns)) {
              desired_columns <- names(x)
            }
            for (col in desired_columns) {
              if (col %in% names(x)) {
                value <- x[[col]]
                if (is.list(value)) {
                  extracted_data[[col]] <- paste(unlist(value), collapse = ",")
                } else {
                  extracted_data[[col]] <- value
                }
              } else {
                extracted_data[[col]] <- NA
              }
            }
            filtered_results <- append(filtered_results, list(extracted_data))
          }, error = function(e) {
            warning(paste("Error processing variant:", e$message))
          })
        }
      }
    }
    
    Sys.sleep(1)  # Add a 1-second delay after each batch request
  }
  
  if (length(filtered_results) > 0) {
    clinvar_df <- do.call(rbind, lapply(filtered_results, as.data.frame, stringsAsFactors = FALSE))
    clinvar_df$uid <- rownames(clinvar_df)
    rownames(clinvar_df) <- NULL
  } else {
    clinvar_df <- NULL
  }
  
  return(clinvar_df)
}

# Example usage:
genes <- c("PTEN")  # You can add more gene names to this vector if needed
desired_cols <- c("name", "clinical_significance", "variation_type", "rsid", "gene", "chromosome", "start", "stop", "obj_type", "molecular_consequence_list", "hgvs", "protein_change")

all_clinvar_data <- data.frame()  # Initialize an empty dataframe to store the results

for (gene_name in genes) {
  clinvar_data <- get_clinvar_variants_by_gene(gene_name, desired_cols)
  if (!is.null(clinvar_data)) {
    print(paste(gene_name, "ClinVar data retrieved successfully:"))
    print(dim(clinvar_data))  # Print the dimensions of the retrieved data
    print(head(clinvar_data))  # Print the first few rows of the retrieved data
    all_clinvar_data <- rbind(all_clinvar_data, clinvar_data)  # Append to the results dataframe
  } else {
    print(paste("Failed to retrieve ClinVar data for", gene_name))
  }
}

# Optionally, you can check the complete results
print("Final combined data:")
print(dim(all_clinvar_data))  # Print the dimensions of the combined data
print(head(all_clinvar_data))  # Print the first few rows of the combined data






















library(httr)
library(XML)

get_clinvar_gene_data <- function(gene_name, desired_columns = NULL, max_retries = 5, delay_seconds = 5, batch_size_efetch = 20) {
  
  # Function to handle API requests with retries and delays
  make_api_request <- function(url, max_retries = 5, delay_seconds = 5) {
    retries <- 0
    success <- FALSE
    response <- NULL
    
    while (retries < max_retries && !success) {
      response <- GET(url)
      if (status_code(response) == 429) {
        print("Rate limit exceeded. Retrying after delay...")
        Sys.sleep(delay_seconds)
        retries <- retries + 1
      } else if (status_code(response) == 200) {
        success <- TRUE
      } else {
        print(paste("Error in request:", status_code(response), url))
        break
      }
    }
    
    if (!success) {
      print(paste("Failed after", max_retries, "retries:", url))
      return(NULL)
    }
    
    return(response)
  }
  
  # 1. Get Variant IDs using esearch (Corrected search strategy)
  encoded_gene_name <- URLencode(gene_name) # Encode the gene name for URL (still good practice)
  esearch_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=", encoded_gene_name, "&retmax=500&retmode=xml") # Search without [gene]
  esearch_response <- make_api_request(esearch_url)
  
  if (is.null(esearch_response)) return(NULL)
  
  esearch_content <- content(esearch_response, "text")
  esearch_xml <- xmlParse(esearch_content)
  
  count_node <- xpathSApply(esearch_xml, "//eSearchResult/count", xmlValue)
  
  total_records <- 0
  if (length(count_node) > 0) {
    total_records <- as.numeric(count_node)
    if (is.na(total_records)) {
      print(paste("Error: Count is not a valid number for", gene_name))
      return(NULL)
    }
  } else {
    print(paste("No count information found for", gene_name, ". Assuming 0 records."))
    return(NULL)
  }
  
  if (total_records == 0) {
    print(paste("No ClinVar records found for gene:", gene_name))
    return(NULL)
  }
  
  print(paste("Total ClinVar records found (before filtering) for", gene_name, ":", total_records))
  
  
  all_ids <- c()
  retmax <- 500
  
  for (retstart in seq(0, total_records, by = retmax)) { # REMOVE MIN in production
    esearch_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=", encoded_gene_name, "&retmax=", retmax, "&retstart=", retstart, "&retmode=xml") # Search without [gene]
    esearch_response <- make_api_request(esearch_url)
    if (is.null(esearch_response)) next
    
    esearch_content <- content(esearch_response, "text")
    esearch_xml <- xmlParse(esearch_content)
    ids <- xpathSApply(esearch_xml, "//IdList/Id", xmlValue)
    all_ids <- c(all_ids, ids)
    Sys.sleep(delay_seconds)
  }
  
  if (length(all_ids) == 0) {
    print(paste("No ClinVar IDs found for gene:", gene_name))
    return(NULL)
  }
  
  # 2. Fetch Variant Details using efetch and FILTER by gene symbol
  all_data <- list()
  for (batch_start in seq(1, length(all_ids), by = batch_size_efetch)) {
    batch_ids <- all_ids[batch_start:min(batch_start + batch_size_efetch - 1, length(all_ids))]
    efetch_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=", paste(batch_ids, collapse = ","), "&from_esearch=true")
    
    efetch_response <- make_api_request(efetch_url)
    if (is.null(efetch_response)) next
    
    efetch_content <- content(efetch_response, "text")
    tryCatch({
      efetch_xml <- xmlParse(efetch_content)
      for (record in getNodeSet(efetch_xml, "//ClinVarVariation")) {
        record_data <- list()
        
        # Extract gene symbol and FILTER
        record_gene_symbol <- xmlValue(record[["Gene"]][["Symbol"]])
        if (record_gene_symbol == gene_name) { # Exact match for gene name
          record_data[["variation_id"]] <- xmlValue(record[["VariationID"]])
          record_data[["clinical_significance"]] <- xmlValue(record[["ClinicalSignificance"]][["Description"]])
          record_data[["gene_symbol"]] <- record_gene_symbol # Store the extracted gene symbol
          record_data[["chromosome"]] <- xmlValue(record[["Location"]][["Chromosome"]])
          record_data[["start"]] <- xmlValue(record[["Location"]][["Start"]])
          record_data[["stop"]] <- xmlValue(record[["Location"]][["Stop"]])
          record_data[["hgvs"]] <- xmlValue(record[["HGVS"]])
          record_data[["protein_change"]] <- xmlValue(record[["ProteinChange"]])
          record_data[["rsid"]] <- xmlValue(record[["RS"]])  # Check if RS exists before extraction
          record_data[["molecular_consequence"]] <- paste(xpathSApply(record, "//MolecularConsequence/Name", xmlValue), collapse=",") # Get all molecular consequences
          
          # ... Add other elements you need ...
          all_data <- append(all_data, list(record_data))
        }
      }
    }, error = function(e) {
      print(paste("Error parsing XML:", e$message))
      print(efetch_content) # Print the problematic XML for debugging
    })
    Sys.sleep(delay_seconds)
  }
  
  if (length(all_data) > 0) {
    clinvar_df <- do.call(rbind, lapply(all_data, as.data.frame, stringsAsFactors = FALSE))
    clinvar_df$uid <- rownames(clinvar_df)
    rownames(clinvar_df) <- NULL
  } else {
    clinvar_df <- NULL
  }
  
  return(clinvar_df)
}

# Example usage:
genes <- c("PTEN") # Example including gene with space
all_gene_data <- data.frame()

for (gene in genes) {
  clinvar_data <- get_clinvar_gene_data(gene, batch_size_efetch = 50) # Adjust batch size as needed
  if (!is.null(clinvar_data)) {
    print(paste(gene, "ClinVar data retrieved successfully:"))
    print(dim(clinvar_data))
    #print(head(clinvar_data)) # Print a few rows for inspection
    all_gene_data <- rbind(all_gene_data, clinvar_data)
  } else {
    print(paste("Failed to retrieve ClinVar data for", gene))
  }
}

print(dim(all_gene_data))
#print(head(all_gene_data)) # Inspect the combined data