### Load libraries 

library(dplyr)
library(ggplot2)
library(pROC)
library(GenomicRanges)
library(Rsamtools)
library(VariantAnnotation)

### Define file paths 

define_regions_and_files <- function() {
  list(
    rng_pten = GRanges(seqnames = "chr10", ranges = IRanges(start = 87863114, end = 87971941, names = "PTEN")),
    rng_CADD = GRanges(seqnames = "10", ranges = IRanges(start = 87863114, end = 87971941, names = "PTEN")),
    gnomad_vcf_chr10 = "F:/Capstone/Resources/gnomAD/gnomad.joint.v4.1.sites.chr10.vcf.bgz",
    clinvar_path = "F:/Capstone/Resources/ClinVar/clinvar.vcf.gz",
    cadd_path = "F:/Capstone/Resources/CADD/v1.7/whole_genome_SNVs.tsv.gz"
  )
} 


### Gnomad data extraction 

extract_gnomad_data <- function(vcf_path, region) {
  tabix <- TabixFile(vcf_path)
  vcf <- readVcf(tabix, "hg38", param = region)
  fixed <- as.data.frame(rowRanges(vcf))
  fixed$FAFmax_faf95_max_joint <- unlist(info(vcf)$fafmax_faf95_max_joint)
  fixed$alt_values <- sapply(fixed$ALT, as.character)
  fixed$merge <- paste(fixed$start, fixed$REF, fixed$alt_values, sep = " ")
  fixed <- fixed %>% filter(FILTER == "PASS")
  return(fixed)
}

### Extract for PTEN 

extract_clinvar_pten_variants <- function(clinvar_path) {
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

### Clinvar/Gnomad Merge

merge_clinvar_gnomad <- function(clinvar_df, gnomad_df) {
  merged <- merge(clinvar_df, gnomad_df, by = "merge", all = TRUE)
  return(merged)
}

### CADD Data extraction 

extract_cadd_data <- function(cadd_path, region) {
  tabix <- TabixFile(cadd_path)
  raw_data <- scanTabix(tabix, param = region)
  df <- read.table(text = unlist(raw_data), header = FALSE, sep = "\t")
  colnames(df) <- c("CHROM", "start", "REF", "ALT", "RawScore", "PHRED")
  df$merge <- paste(df$start, df$REF, df$ALT, sep = " ")
  return(df)
}

### CADD Merge 

merge_all_annotations <- function(merged_df, cadd_df) {
  combined <- merge(merged_df, cadd_df, by = "merge", all.x = TRUE)
  return(combined)
}

### Pathogenicity assignment 

assign_pathogenic_labels <- function(df, faf_cutoff = 10^-4) {
  df <- df %>%
    mutate(
      CLNSIG = as.character(CLNSIG),
      CLNSIG = ifelse(CLNSIG %in% c("NULL", "character(0)"), NA, CLNSIG),
      pathogenic = case_when(
        CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic") ~ 1,
        CLNSIG %in% c("Benign", "Likely_benign", "Benign/Likely_benign") ~ -1,
        is.na(CLNSIG) & is.na(FAFmax_faf95_max_joint) ~ 0,
        is.na(CLNSIG) & FAFmax_faf95_max_joint > faf_cutoff ~ -1,
        TRUE ~ 0
      )
    )
  return(df)
}

### Weighted logistic regression 

prob_imputation <- function(filtered_data, prop_path) {
  total_obs <- sum(!is.na(filtered_data$pathogenic))
  n_case <- sum(filtered_data$pathogenic == 1, na.rm = TRUE)
  n_cont <- sum(filtered_data$pathogenic == -1, na.rm = TRUE)
  
  w_case <- total_obs * prop_path / n_case
  w_cntl <- total_obs * (1 - prop_path) / n_cont
  
  d_observed <- filtered_data %>%
    filter(pathogenic != 0) %>%
    mutate(w = ifelse(pathogenic == 1, w_case, w_cntl))
  
  d_unc <- filtered_data %>%
    filter(pathogenic == 0)
  
  model_estimates <- list()
  model_estimates_weight <- list()
  plot_data <- c()
  
  for (annot in c("PHRED")) {
    d_obs <- d_observed %>%
      select(merge, all_of(annot), pathogenic, w) %>%
      mutate(pathogenic = as.numeric(ifelse(pathogenic == -1, 0, 1)))
    
    formula <- as.formula(paste("pathogenic ~", annot))
    
    model <- glm(formula, data = d_obs, family = "binomial", weights = w)
    model_noweight <- glm(formula, data = d_obs, family = "binomial")
    
    d_unc <- filtered_data %>%
      filter(pathogenic == 0) %>%
      select(merge, pathogenic, all_of(annot)) %>%
      mutate(prob_disease = predict(model, newdata = ., type = "response")) %>%
      select(merge, prob_disease)
    
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
        plot_data <- d_imp %>% select(pathogenic, all_of(annot))
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
    model_estimates_weight <- list(params = summary(model)$coef[,1],
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

### Threshold Table 

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


### RUN FULL PIPELINE 

run_pipeline <- function() {
  paths <- define_regions_and_files()
  
  gnomad <- extract_gnomad_data(paths$gnomad_vcf_chr10, paths$rng_pten)
  clinvar <- extract_clinvar_pten_variants(paths$clinvar_path)
  merged <- merge_clinvar_gnomad(clinvar, gnomad)
  cadd <- extract_cadd_data(paths$cadd_path, paths$rng_CADD)
  combined <- merge_all_annotations(merged, cadd)
  labeled <- assign_pathogenic_labels(combined)
  filtered_data <- labeled %>% filter(!is.na(PHRED))
  
  prop_path <- 0.04 # argument into pipeline (set this as default but have user be able to input a value if they want)
  tmp <- prob_imputation(filtered_data, prop_path)
  threshold_table <- calculate_thresholds(tmp$model_estimates, filtered_data, prop_path)
  
  print(threshold_table)
} 

run_pipeline()