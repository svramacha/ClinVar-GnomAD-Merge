# Load necessary libraries
library(GenomicRanges)
library(httr)
library(jsonlite)
library(Rsamtools)
library(VariantAnnotation)

# Function to combine gnomAD, ClinVar, and CADD data for PTEN gene on chr10
combine_pten_data <- function(gnomad_url, clinvar_vcf_path, cadd_file_path, chrom = "chr10", start = 87863114, end = 87971941) {
  
  # Define the region of interest for PTEN (chromosome 10)
  rng_pten <- GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end, names = "PTEN"))
  
  # ---- 1. Query gnomAD API for PTEN gene data on chromosome 10 ----
  query <- list(
    "gene" = "PTEN",
    "chrom" = chrom,
    "start" = as.character(start),
    "end" = as.character(end)
  )
  
  # Send the API request to gnomAD
  response <- GET(gnomad_url, query = query)
  
  # Parse the response
  response_data <- content(response, "text")
  gnomad_data <- fromJSON(response_data)
  
  # Extract relevant information from the gnomAD response
  variants_data <- gnomad_data$data$variants
  
  # Create a dataframe for gnomAD variants
  gnomad_df <- data.frame(
    CHROM = chrom,
    POS = unlist(lapply(variants_data, function(x) x$start)),
    REF = unlist(lapply(variants_data, function(x) x$reference)),
    ALT = unlist(lapply(variants_data, function(x) x$alt)),
    AF = unlist(lapply(variants_data, function(x) x$allele_frequency)),
    stringsAsFactors = FALSE
  )
  
  # ---- 2. Extract ClinVar Data for PTEN ----
  vcf_clinvar <- readVcf(clinvar_vcf_path, "hg38")
  
  # Filter for PTEN variants in ClinVar
  geneinfo_data <- info(vcf_clinvar)$GENEINFO
  pten_variants <- grepl("PTEN", geneinfo_data, ignore.case = TRUE)
  vcf_pten <- vcf_clinvar[pten_variants, ]
  
  # Extract information from ClinVar VCF
  info_pten <- info(vcf_pten)
  pten_info <- as.data.frame(info_pten)
  
  gr_pten <- rowRanges(vcf_pten)
  pten_gr <- as.data.frame(gr_pten)
  
  pten_df <- cbind(pten_info, pten_gr)
  
  # ---- 3. Query CADD data for PTEN region ----
  # Open the CADD file with TabixFile
  tabix_file <- TabixFile(cadd_file_path)
  
  # Query the CADD file for the PTEN region on chromosome 10
  cadd_data <- scanTabix(tabix_file, param = rng_pten)
  
  # Convert the extracted data into a dataframe
  cadd_df <- read.table(text = unlist(cadd_data), header = FALSE, sep = "\t")
  
  # Change the column names so that they match for merging
  colnames(cadd_df) <- c("CHROM", "start", "REF", "ALT", "RawScore", "PHRED")
  
  # Create merge column for CADD
  cadd_df$merge <- paste(cadd_df$start, cadd_df$REF, cadd_df$ALT, sep = " ")
  
  # ---- 4. Merge gnomAD, ClinVar, and CADD data ----
  # Prepare ClinVar data
  pten_df$alt_values <- sapply(pten_df$ALT, function(x) as.character(x))
  pten_df$merge <- paste(pten_df$start, pten_df$REF, pten_df$alt_values, sep = " ")
  
  # Merge gnomAD and ClinVar data (outer join)
  combined_data <- merge(pten_df, gnomad_df, by = "merge", all = TRUE)
  
  # Merge with CADD data (left join to preserve all ClinVar + gnomAD data)
  final_combined_data <- merge(combined_data, cadd_df, by = "merge", all.x = TRUE)
  
  # Return the final combined data
  return(final_combined_data)
}

# Example Usage:
gnomad_url <- "https://gnomad.broadinstitute.org/api/variant/"
clinvar_vcf_path <- "F:/Capstone/Resources/ClinVar/clinvar.vcf.gz"
cadd_file_path <- "F:/Capstone/Resources/CADD/v1.7/whole_genome_SNVs.tsv.gz"

# Get the combined data for PTEN gene on chr10
pten_combined_data <- combine_pten_data(gnomad_url, clinvar_vcf_path, cadd_file_path)

# Print the resulting dataframe
print(pten_combined_data)
