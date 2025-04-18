### Still need to double check o doing the multiple inputs ? or have a user input available - automating the regions

# ------------------------------
# CHUNK 1: Load Libraries
# ------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
library(GenomicRanges)
library(Rsamtools)
library(VariantAnnotation)

# ------------------------------
# CHUNK 2: Extract gnomAD for PTEN
# Function extracts variant data from the gnomAD VCF file for the PTEN gene region.
# Constructs a GRanges object to define the genomic interval and returns a cleaned dataframe
# with key fields including allele frequency (FAF) and a unique merge identifier.
# ------------------------------
extract_gnomad_pten <- function(vcf_path) {
  region <- GRanges(seqnames = "chr10", ranges = IRanges(start = 87863114, end = 87971941, names = "PTEN"))
  tabix <- TabixFile(vcf_path)
  vcf <- readVcf(tabix, "hg38", param = region)
  
  df <- as.data.frame(rowRanges(vcf))
  df$FAFmax_faf95_max_joint <- unlist(info(vcf)$fafmax_faf95_max_joint)
  df$alt_values <- sapply(df$ALT, as.character)
  df$merge <- paste(df$start, df$REF, df$alt_values, sep = " ")
  df <- df %>% filter(FILTER == "PASS")
  
  return(df)
}

# ------------------------------
# CHUNK 3: Extract ClinVar for PTEN
# Function reads ClinVar VCF, extracts all PTEN-annotated variants using the GENEINFO field,
# and returns a dataframe with INFO fields and genomic ranges for merging.
# ------------------------------
extract_clinvar_pten <- function(clinvar_path) {
  vcf <- readVcf(clinvar_path, "hg38")
  geneinfo <- info(vcf)$GENEINFO
  pten_variants <- grepl("PTEN", geneinfo, ignore.case = TRUE)
  vcf_pten <- vcf[pten_variants, ]
  
  info_df <- as.data.frame(info(vcf_pten))
  gr_df <- as.data.frame(rowRanges(vcf_pten))
  df <- cbind(info_df, gr_df)
  df$alt_values <- sapply(df$ALT, as.character)
  df$merge <- paste(df$start, df$REF, df$alt_values, sep = " ")
  
  return(df)
}

# ------------------------------
# CHUNK 4: Extract CADD for PTEN
# Function uses tabix to scan the CADD file and extract annotation scores for PTEN region.
# CADD scores like PHRED and RawScore help estimate deleteriousness of variants.
# ------------------------------
extract_cadd_pten <- function(cadd_path) {
  region <- GRanges(seqnames = "10", ranges = IRanges(start = 87863114, end = 87971941, names = "PTEN"))
  tabix <- TabixFile(cadd_path)
  raw_data <- scanTabix(tabix, param = region)
  
  df <- read.table(text = unlist(raw_data), header = FALSE, sep = "\t")
  colnames(df) <- c("CHROM", "start", "REF", "ALT", "RawScore", "PHRED")
  df$merge <- paste(df$start, df$REF, df$ALT, sep = " ")
  
  return(df)
}

# ------------------------------
# CHUNK 5: Merge All Annotations
# Merges ClinVar, gnomAD, and CADD data using the merge key (start, REF, ALT).
# This prepares a unified dataset for pathogenicity classification and modeling.
# ------------------------------
merge_all_annotations <- function(clinvar_df, gnomad_df, cadd_df) {
  merged <- merge(clinvar_df, gnomad_df, by = "merge", all = TRUE)
  merged_final <- merge(merged, cadd_df, by = "merge", all.x = TRUE)
  return(merged_final)
}

# ------------------------------
# CHUNK 6–8: Modeling and Thresholds
# This function performs the full pipeline of logistic regression modeling using PHRED scores.
# It includes probabilistic imputation for uncertain variants and calculates ACMG score thresholds.
# Inputs:
#   df = merged dataframe from ClinVar + gnomAD + CADD
#   prop_path = estimated proportion of pathogenic variants
#   faf_cutoff = threshold for gnomAD FAF to define benign variants
#   threshold_output = output file name for threshold table
# Output:
#   A list with filtered data, model results, and ACMG threshold table
# ------------------------------

run_analysis <- function(df, prop_path = 0.04, faf_cutoff = 10^-4, threshold_output = "acmg_table_impute.txt") {
  # STEP 1: Assign pathogenicity from ClinVar and gnomAD
  # 1 = Pathogenic, -1 = Benign, 0 = Uncertain/VUS
  df <- df %>%
    mutate(CLNSIG = as.character(CLNSIG),
           CLNSIG = ifelse(CLNSIG %in% c("NULL", "character(0)"), NA, CLNSIG)) %>%
    mutate(pathogenic = case_when(
      CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic") ~ 1,
      CLNSIG %in% c("Benign", "Likely_benign", "Benign/Likely_benign") ~ -1,
      is.na(CLNSIG) & is.na(FAFmax_faf95_max_joint) ~ 0,
      is.na(CLNSIG) & FAFmax_faf95_max_joint > faf_cutoff ~ -1,
      TRUE ~ 0
    ))
  
  # Filter for records with PHRED scores
  filtered_data <- df %>% filter(!is.na(PHRED))
  
  # These include logistic regression with weights, 25x imputation, and uniroot thresholding
  prob_imputation <- function(filtered_data, prop_path) {
    total_obs <- sum(!is.na(filtered_data$pathogenic))
    n_case <- sum(filtered_data$pathogenic == 1, na.rm = TRUE)
    n_cont <- sum(filtered_data$pathogenic == -1, na.rm = TRUE)
    
    w_case <- total_obs * prop_path / n_case
    w_cntl <- total_obs * (1 - prop_path) / n_cont
    
    d_observed <- filtered_data %>%
      filter(pathogenic != 0) %>%
      mutate(w = ifelse(pathogenic == 1, w_case, w_cntl))
    
    model_estimates <- list()
    model_estimates_weight <- list()
    plot_data <- c()
    
    for (annot in c("PHRED")) {
      d_obs <- d_observed %>%
        dplyr::select(merge, all_of(annot), pathogenic, w) %>%
        mutate(pathogenic = as.numeric(ifelse(pathogenic == -1, 0, 1)))
      
      formula <- as.formula(paste("pathogenic ~", annot))
      
      model <- glm(formula, data = d_obs, family = "binomial", weights = w)
      
      d_unc <- filtered_data %>%
        filter(pathogenic == 0) %>%
        dplyr::select(merge, pathogenic, all_of(annot)) %>%
        mutate(prob_disease = predict(model, newdata = ., type = "response")) %>%
        dplyr::select(merge, prob_disease)
      
      d_all <- filtered_data %>%
        mutate(weight = 1) %>%
        left_join(d_unc, by = "merge") %>%
        mutate(weight = ifelse(pathogenic == 0, abs(prob_disease - 0.5) * 2, weight))
      
      n_imputations <- 25
      models <- vector("list", n_imputations)
      
      for (i in 1:n_imputations) {
        d_imp <- d_all %>%
          rowwise() %>%
          mutate(pathogenic = ifelse(!is.na(prob_disease), rbinom(1, 1, prob_disease), pathogenic),
                 pathogenic = ifelse(pathogenic == -1, 0, pathogenic)) %>%
          ungroup()
        
        models[[i]] <- glm(formula, family = binomial, data = d_imp, weights = weight)
        
        if (i == 1) {
          plot_data <- d_imp %>% dplyr::select(pathogenic, all_of(annot))
        }
      }
      
      vcov_within <- Reduce("+", lapply(models, vcov)) / n_imputations
      beta_bar <- rowMeans(sapply(models, coef))
      vcov_between <- Reduce("+", lapply(models, function(m) {
        diff <- coef(m) - beta_bar
        tcrossprod(diff)
      })) / (n_imputations - 1)
      
      vcov_total <- vcov_within + (1 + 1 / n_imputations) * vcov_between
      t_stat <- beta_bar / sqrt(diag(vcov_total))
      p_value <- 2 * pt(-abs(t_stat), df = n_imputations - 1)
      
      model_estimates[[annot]] <- list(params = beta_bar, cov = vcov_total, t = t_stat, p = p_value)
      model_estimates_weight[[annot]] <- list(params = summary(model)$coef[,1],
                                              cov = vcov(model),
                                              t = summary(model)$coef[,3],
                                              p = summary(model)$coef[,4])
    }
    
    return(list(
      model_estimates = model_estimates,
      model_estimates_weight = model_estimates_weight,
      aucs = NULL,
      roc_data = NULL,
      plot_data = plot_data
    ))
  }
  
  calculate_thresholds <- function(est_perm, filtered_data, prop_path, file = "acmg_table_impute.txt") {
    rslts <- c()
    cutoffs <- c(2.406, 5.790, 33.53, 1124)
    cutlabel <- c("supporting", "moderate", "strong", "very_strong")
    odds_path <- prop_path / (1 - prop_path)
    
    min_func <- function(x, beta, v, cutoff, odds_path, var_type) {
      var_adj <- ifelse(var_type == "path", 1, -1)
      se_log_odds <- sqrt(t(c(1, x)) %*% v %*% c(1, x))
      z <- qnorm(.95)
      logodds_bound <- (beta[1] + beta[2] * x) - var_adj * z * se_log_odds
      odds_bound <- exp(var_adj * logodds_bound)
      odds_path_adj <- odds_path^var_adj
      return(cutoff - odds_bound / odds_path_adj)
    }
    
    for (annot in c("PHRED")) {
      betas <- est_perm[[annot]]$params
      v <- est_perm[[annot]]$cov
      int <- betas[1]
      logor <- betas[2]
      p <- est_perm[[annot]]$p[2]
      se_int <- sqrt(v[1, 1])
      se_logor <- sqrt(v[2, 2])
      rng <- range(filtered_data[[annot]], na.rm = TRUE)
      lower <- min(rng)
      upper <- max(rng)
      
      for (j in seq_along(cutoffs)) {
        cutoff <- cutoffs[j]
        
        result_path <- tryCatch({
          uniroot(min_func, lower = lower, upper = upper, beta = betas, v = v,
                  cutoff = cutoff, odds_path = odds_path, var_type = "path")
        }, error = function(e) list(root = NA))
        
        result_benign <- tryCatch({
          uniroot(min_func, lower = lower, upper = upper, beta = betas, v = v,
                  cutoff = cutoff, odds_path = odds_path, var_type = "benign")
        }, error = function(e) list(root = NA))
        
        rslts <- rbind(rslts,
                       c(annot, int, logor, se_int, se_logor, p, "pathogenic", cutlabel[j], result_path$root, rng),
                       c(annot, int, logor, se_int, se_logor, p, "benign", cutlabel[j], result_benign$root, rng))
      }
    }
    
    rslts <- as.data.frame(rslts, stringsAsFactors = FALSE)
    names(rslts) <- c("annot", "int", "logor", "se_int", "se_logor", "pval", "var_type",
                      "acmg_cat", "annot_value", "min_obs", "max_obs")
    
    tbl <- rslts %>%
      tidyr::pivot_wider(names_from = c("var_type", "acmg_cat"),
                         values_from = "annot_value",
                         names_glue = "{var_type}_{acmg_cat}") %>%
      dplyr::select(annot, int, logor, se_int, se_logor, pval,
                    starts_with("pathogenic"), starts_with("benign"),
                    min_obs, max_obs)
    
    write.table(file = file, tbl, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    cat("Wrote ACMG threshold table to", file, "\n")
    
    return(tbl)
  }
  
  tmp <- prob_imputation(filtered_data, prop_path)
  threshold_table <- calculate_thresholds(tmp$model_estimates, filtered_data, prop_path, file = threshold_output)
  
  return(list(
    data = filtered_data,
    models = tmp,
    threshold_table = threshold_table
  ))
}


# ------------------------------
# CHUNK 9: Run Entire Pipeline for PTEN
# ------------------------------

# Define file paths to VCF and annotation sources
gnomad_path <- "F:/Capstone/Resources/gnomAD/gnomad.joint.v4.1.sites.chr10.vcf.bgz"
clinvar_path <- "F:/Capstone/Resources/ClinVar/clinvar.vcf.gz"
cadd_path <- "F:/Capstone/Resources/CADD/v1.7/whole_genome_SNVs.tsv.gz"

# Extract variant data for PTEN from each source
gnomad_df <- extract_gnomad_pten(gnomad_path)
clinvar_df <- extract_clinvar_pten(clinvar_path)
cadd_df <- extract_cadd_pten(cadd_path)

# Merge all sources into one combined dataframe
final_df <- merge_all_annotations(clinvar_df, gnomad_df, cadd_df)

# Run model + thresholds
results <- run_analysis(final_df)
print(results$threshold_table)
