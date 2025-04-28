library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
library(GenomicRanges)
library(Rsamtools)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)

# ------------------------
# Extract Gnomad via API
#-------------------------

get_variants_by_gene <- function(gene_name) {
  library(httr)
  library(jsonlite)
  library(dplyr)
  
  # Step 1: Get gene coordinates using gnomAD API
  gene_query <- '
  query GeneInfo($gene_name: String!) {
    gene(gene_symbol: $gene_name, reference_genome: GRCh38) {
      chrom
      start
      stop
    }
  }'
  
  variables <- list(gene_name = gene_name)
  
  gene_body <- toJSON(list(query = gene_query, variables = variables), auto_unbox = TRUE)
  url <- "https://gnomad.broadinstitute.org/api"
  
  gene_response <- POST(url, body = gene_body, encode = "raw", content_type("application/json"))
  
  if (status_code(gene_response) != 200) {
    stop("Error fetching gene info.")
  }
  
  gene_data <- fromJSON(content(gene_response, "text"))
  region <- gene_data$data$gene
  
  if (is.null(region)) {
    stop(paste("Gene", gene_name, "not found in gnomAD."))
  }
  
  # Step 2: Query variants in that region
  variant_query <- '
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
  }'
  
  variant_vars <- list(
    chrom = region$chrom,
    start = region$start,
    stop = region$stop
  )
  
  variant_body <- toJSON(list(query = variant_query, variables = variant_vars), auto_unbox = TRUE)
  
  variant_response <- POST(url, 
                           body = variant_body, 
                           encode = "raw", 
                           content_type("application/json"),
                           timeout(60))
  
  if (status_code(variant_response) != 200) {
    stop("Error fetching variant data.")
  }
  
  variant_data <- fromJSON(content(variant_response, "text"))
  variants <- variant_data$data$region$variants
  
  if (length(variants) == 0) {
    message("No variants found for ", gene_name)
    return(data.frame())
  }
  
  # Clean up the results
  variants_df <- as.data.frame(variants)
  
  variants_df$flags <- sapply(variants_df$flags, paste, collapse = ",")
  variants_df$filters <- sapply(variants_df$joint$filters, paste, collapse = ",")
  variants_df$ac <- variants_df$joint$ac
  variants_df$an <- variants_df$joint$an
  variants_df$faf95_max <- variants_df$joint$fafmax$faf95_max
  variants_df$faf95_max_gen_anc <- variants_df$joint$fafmax$faf95_max_gen_anc
  variants_df$joint <- NULL
  
  variants_df$chrom <- region$chrom
  variants_df$start <- region$start
  variants_df$stop <- region$stop
  
  return(variants_df)
}

# Example usage
my_gene <- "PTEN"
variants <- get_variants_by_gene(my_gene)
head(variants)

any(duplicated(variants$variant_id)) #no duplicating variant_id's





# ------------------------
# Extract Clinvar via API
#-------------------------

get_clinvar_variants <- function(gene_name) {
  library(httr)
  library(jsonlite)
  library(dplyr)
  
  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  
  # Step 1: Get total number of results
  count_url <- paste0(
    base_url,
    "esearch.fcgi?db=clinvar&term=", gene_name, "[gene]&retmode=json"
  )
  
  response <- GET(count_url)
  ids_data <- fromJSON(content(response, "text"))
  total_records <- as.numeric(ids_data$esearchresult$count)
  
  if (total_records == 0) {
    message("No ClinVar records found for ", gene_name)
    return(data.frame())
  }
  
  # Step 2: Fetch all UIDs
  esearch_url <- paste0(
    base_url,
    "esearch.fcgi?db=clinvar&term=", gene_name, "[gene]",
    "&retmax=", total_records, "&retmode=json"
  )
  
  response <- GET(esearch_url)
  ids_data <- fromJSON(content(response, "text"))
  uids <- ids_data$esearchresult$idlist
  
  # Step 3: Get variant summaries
  chunk_size <- 200  # NCBI recommends batching requests
  result_list <- list()
  
  for (i in seq(1, length(uids), by = chunk_size)) {
    chunk_ids <- uids[i:min(i + chunk_size - 1, length(uids))]
    ids_str <- paste(chunk_ids, collapse = ",")
    
    esummary_url <- paste0(
      base_url,
      "esummary.fcgi?db=clinvar&id=", ids_str, "&retmode=json"
    )
    
    summary_response <- GET(esummary_url)
    summary_data <- fromJSON(content(summary_response, "text"))
    
    chunk_results <- summary_data$result[chunk_ids]
    
    list_of_records <- lapply(chunk_results, function(x) {
      flat <- unlist(x, recursive = TRUE)
      as.data.frame(as.list(flat), stringsAsFactors = FALSE)
    })
    
    result_list <- c(result_list, list_of_records)
  }
  
  df <- bind_rows(result_list)
  return(df)
}


# Example usage 
clinvar_df <- get_clinvar_variants("PTEN")  # Use the relevant gene



# Check unique values in the 'obj_type' column
unique_obj_types <- unique(clinvar_df$obj_type)

# Count the number of unique values
num_unique_obj_types <- length(unique_obj_types)

# Print the result
print(num_unique_obj_types) #about 2600 SNV variants


# # Verify the count of SNVs specifically
# gene_name <- "PTEN"
# 
# snv_count_url <- paste0(
#   base_url,
#   "esearch.fcgi?db=clinvar&term=SNV[variation_type]+AND+", gene_name, "[gene]&retmode=json"
# )
# 
# response <- GET(snv_count_url)
# snv_data <- fromJSON(content(response, "text"))
# total_snv_records <- as.numeric(snv_data$esearchresult$count)
# 
# print(total_snv_records)




