#### CHUNK 1: Libraries 

library(dplyr)
library(ggplot2)
library(pROC)
library(GenomicRanges)
library(Rsamtools)
library(VariantAnnotation)


#### CHUNK 2: Define regions and paths 

# Define region of interest
rng <- GRanges(seqnames="chr10", ranges=IRanges(
  start=c(89765432), 
  end=c(89774321),
  names=c("gene_123456")))

# Define region of PTEN: 
rng_pten <- GRanges(seqnames = "chr10", ranges = IRanges(
  start = 87863114,
  end = 87971941,
  names = "PTEN"
))

# Define region for CADD: 
rng_CADD <- GRanges(seqnames = "10", ranges = IRanges(
  start = 87863114,
  end = 87971941,
  names = "PTEN"
))

# Path to gnomAD VCF for chromosome 10
gnomad_vcf_chr10 <- "F:/Capstone/Resources/gnomAD/gnomad.joint.v4.1.sites.chr10.vcf.bgz"

#### CHUNK 3: Extract Gnomad 

tab_gnomad <- TabixFile(gnomad_vcf_chr10)
vcf_rng_gnomad_pten <- readVcf(tab_gnomad, "hg38", param=rng_pten)

gnomad_fixed <- as.data.frame(rowRanges(vcf_rng_gnomad_pten))

info <- info(vcf_rng_gnomad_pten)
faf <- unlist(info$fafmax_faf95_max_joint)

gnomad_fixed$FAFmax_faf95_max_joint <- faf
gnomad_fixed$alt_values <- sapply(gnomad_fixed$ALT, function(x) as.character(x))
gnomad_fixed$merge <- paste(gnomad_fixed$start, gnomad_fixed$REF, gnomad_fixed$alt_values, sep = " ")
gnomad_fixed <- gnomad_fixed %>% filter(FILTER == "PASS")

#### CHUNK 4: Extract for PTEN 

clinvar_vcf <- 'F:/Capstone/Resources/ClinVar/clinvar.vcf.gz'
vcf_clinvar <- readVcf(clinvar_vcf, "hg38")

geneinfo_data <- info(vcf_clinvar)$GENEINFO
pten_variants <- grepl("PTEN", geneinfo_data, ignore.case = TRUE)
vcf_pten <- vcf_clinvar[pten_variants, ]

info_pten <- info(vcf_pten)
pten_info <- as.data.frame(info_pten)
gr_pten <- rowRanges(vcf_pten)
pten_gr <- as.data.frame(gr_pten)
pten_df <- cbind(pten_info, pten_gr)
pten_df$alt_values <- sapply(pten_df$ALT, function(x) as.character(x))
pten_df$merge <- paste(pten_df$start, pten_df$REF, pten_df$alt_values, sep = " ")

#### CHUNK 5: Merge Clinvar and Gnomad 

combined_outer_join_pten <- merge(pten_df, gnomad_fixed, by = "merge", all = TRUE)

#### CHUNK 6: Extract CADD

cadd_file <- "F:/Capstone/Resources/CADD/v1.7/whole_genome_SNVs.tsv.gz"
tabix_file <- TabixFile(cadd_file)
cadd_data <- scanTabix(tabix_file, param = rng_CADD)

cadd_df <- read.table(text = unlist(cadd_data), header = FALSE, sep = "\t")
colnames(cadd_df) <- c("CHROM", "start", "REF", "ALT", "RawScore", "PHRED")
cadd_df$merge <- paste(cadd_df$start, cadd_df$REF, cadd_df$ALT, sep = " ")

#### CHUNK 7: Final merge of DB 

final_combined_data_2 <- merge(
  combined_outer_join_pten, 
  cadd_df, 
  by = "merge", 
  all.x = TRUE
)

#### CHUNK 8: Assign pathogenicty 

faf_cutoff <- 10^-4

final_combined_data_2 <- final_combined_data_2 %>% mutate(CLNSIG = as.character(CLNSIG), 
                                                          CLNSIG = ifelse(CLNSIG %in% c("NULL", "character(0)"), NA, CLNSIG))

final_combined_data_2 <- final_combined_data_2 %>%
  mutate(
    pathogenic = ifelse(
      CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"), 1,
      ifelse(
        CLNSIG %in% c("Benign", "Likely_benign", "Benign/Likely_benign"), -1,
        ifelse(
          is.na(CLNSIG) & is.na(FAFmax_faf95_max_joint), 0,
          ifelse(
            is.na(CLNSIG) & FAFmax_faf95_max_joint > faf_cutoff, -1, 0
          )
        )
      )
    )
  )

#### CHUNK 9: double checking info - do we need this? 

table(final_combined_data_2$pathogenic)
sum(is.na(final_combined_data_2$pathogenic))

#### CHUNK 10: First Log Model 

# Filter to rows where PHRED is not NA and pathogenic is not 0 (remove VUS)
filtered_data <- final_combined_data_2 %>% filter(!is.na(PHRED))
path_benign <- filtered_data %>%
  filter(pathogenic != 0) %>%
  mutate(pathogenic = factor(pathogenic, levels = c("-1", "1")))

# Run logistic regression
logistic_model <- glm(pathogenic ~ PHRED, data = path_benign, family = binomial())

# View summary
summary(logistic_model)

# Optional: visualize
boxplot(PHRED ~ pathogenic, data = path_benign, main = "PHRED by Pathogenicity")

#### CHUNK 11: Weighted Log Model with Imputations (already a fxn)

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

# Run the function
prop_path_update <- 0.04
tmp <- prob_imputation(filtered_data, prop_path_update)


#### CHUNK 12: ACMG Threshold Tables (already in a fxn)

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

# Run it
threshold_table <- calculate_thresholds(tmp$model_estimates, filtered_data, prop_path_update)
print(threshold_table)

##### LAST CHUNK - LOOK FOR VISUALS 


