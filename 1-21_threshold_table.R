# Building logistic regression model with chromosome 10

library(dplyr)
library(ggplot2)
library(pROC)
library(GenomicRanges)
library(Rsamtools)
library(VariantAnnotation)



# Define region of interest
# Specific for CHR 10
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

# Open the VCF with TabixFile for subsetting
tab_gnomad <- TabixFile(gnomad_vcf_chr10)
vcf_rng_gnomad_pten <- readVcf(tab_gnomad, "hg38", param=rng_pten)

# Extract fixed fields for CHROM, POS, REF, ALT
gnomad_fixed <- as.data.frame(rowRanges(vcf_rng_gnomad_pten))

## Add in allele frequency (the faf95_max_joint) 
info <- info(vcf_rng_gnomad_pten)
faf <- unlist(info$fafmax_faf95_max_joint)

# Extract the FAFmax_faf95_max_joint allele frequency from INFO fields
gnomad_fixed$FAFmax_faf95_max_joint <- faf

# Convert the ALT column to character
gnomad_fixed$alt_values <- sapply(gnomad_fixed$ALT, function(x) as.character(x))

# Get unique values from the new column
unique_alt_values <- unique(gnomad_fixed$alt_values)

# Print the unique values
print(unique_alt_values)

gnomad_fixed$merge <- paste(gnomad_fixed$start, gnomad_fixed$REF, gnomad_fixed$alt_values, sep = " ")

# Keep only FILTER = PASS
gnomad_fixed <- gnomad_fixed %>% 
  filter(FILTER == "PASS")

# Add in Clinvar data 
clinvar_vcf <- 'F:/Capstone/Resources/ClinVar/clinvar.vcf.gz'
vcf_clinvar <- readVcf(clinvar_vcf, "hg38")

geneinfo_data <- info(vcf_clinvar)$GENEINFO

contains_pten <- grep("PTEN", geneinfo_data, ignore.case = TRUE, value = TRUE)

pten_variants <- grepl("PTEN", geneinfo_data, ignore.case = TRUE)

vcf_pten <- vcf_clinvar[pten_variants, ]

info_pten <- info(vcf_pten)
pten_info <- as.data.frame(info_pten)

gr_pten <- rowRanges(vcf_pten)
pten_gr <- as.data.frame(gr_pten)

pten_df <- cbind(pten_info, pten_gr)

### MERGE CLINVAR AND GNOMAD 
# Convert the ALT column to character for ClinVar
pten_df$alt_values <- sapply(pten_df$ALT, function(x) as.character(x))

# Create new column 
pten_df$merge <- paste(pten_df$start, pten_df$REF, pten_df$alt_values, sep = " ")

# Outer join to merge Clivar and Gnomad
combined_outer_join_pten <- merge(pten_df, gnomad_fixed, by = "merge", all = TRUE)

# Adding CADD
cadd_file <- "F:/Capstone/Resources/CADD/v1.7/whole_genome_SNVs.tsv.gz"

# Open the CADD file with TabixFile
tabix_file <- TabixFile(cadd_file)

# Query the CADD file for ch10
cadd_data <- scanTabix(tabix_file, param = rng_CADD)

# Convert the extracted data into a dataframe
cadd_df <- read.table(text = unlist(cadd_data), header = FALSE, sep = "\t")

## Change the column names so that merge is possible 
colnames(cadd_df) <- c("CHROM", "start", "REF", "ALT", "RawScore", "PHRED")
head(cadd_df)

# Create merge column for CADD 
cadd_df$merge <- paste(cadd_df$start, cadd_df$REF, cadd_df$ALT, sep = " ")

# Combined with inner join and CADD 
final_combined_data_2 <- merge(
  combined_outer_join_pten, 
  cadd_df, 
  by = "merge", 
  all.x = TRUE  # Keep all rows from combined_inner_join
)

# Adding in logistic regression 

# Check data summary
table(final_combined_data_2$pathogenic)
str(final_combined_data_2$pathogenic)

# Same assingment as above with faf_cutoff (didn't run)
faf_cutoff <- 10^-4
prop_path <- 0.4

final_combined_data_2 <- final_combined_data_2 %>% mutate(CLNSIG = as.character(CLNSIG), 
                                                          CLNSIG = ifelse(CLNSIG %in% c("NULL", "character(0)"), NA, CLNSIG))

final_combined_data_2 <- final_combined_data_2 %>%
  mutate(
    pathogenic = ifelse(
      CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"), 
      1,  # Assign 1 for pathogenic classifications
      ifelse(
        CLNSIG %in% c("Benign", "Likely_benign", "Benign/Likely_benign"), 
        -1,  # Assign -1 for benign classifications
        ifelse(
          is.na(CLNSIG) & is.na(FAFmax_faf95_max_joint), 
          0,  # Assign 0 if both CLNSIG and FAFmax_faf95_max_joint are NA
          ifelse(
            is.na(CLNSIG) & FAFmax_faf95_max_joint > faf_cutoff, 
            -1,  # Assign -1 if FAFmax_faf95_max_joint exceeds faf_cutoff and CLNSIG is NA
            0  # Default to 0 for all other conditions
          )
        )
      )
    )
  )

table(final_combined_data_2$pathogenic, final_combined_data_2$CLNSIG, useNA = "ifany") #CLNSIG shouldn't be a list
table(final_combined_data_2$pathogenic, final_combined_data_2$FILTER.y)

# TO DO 1/31:
# Assign the one's w/o faf_max or CLNSIG as 0 -- done
# Only keep FILTER = PASS in gnomad_fixed (check if these with > faf_cutoff end up as 0) -- filtered but unsure if they assigned as 0
filter_path <- final_combined_data_2 %>%
  dplyr::filter(FILTER.y == "PASS", FAFmax_faf95_max_joint > faf_cutoff) %>%
  dplyr::select(merge, pathogenic)

table(filter_path$pathogenic)


sum(is.na(final_combined_data_2$pathogenic))

table(final_combined_data_2$pathogenic)
boxplot(PHRED ~ pathogenic, data = path_benign, main = "PHRED by Pathogenicity")

# First regression preformed for benign/pathogenic - have to create a separate df for this regression 
# create lm model data 
# take out anything = 0 
lm_data <- final_combined_data_2 %>%
  filter(pathogenic != 0)

# Build a logistic regression model using CADD PHRED score
# Filter the dataset - filter this before hand in the final combined data part (move this up)
filtered_data <- final_combined_data_2 %>%
  filter(!is.na(PHRED)) 

# Create new dataset (path/benign) filtering out 0's
path_benign <- filtered_data %>%
  filter(pathogenic != 0) %>% mutate(pathogenic=factor(pathogenic, levels = c("-1","1"))) 

logistic_model <- glm(pathogenic ~ PHRED, data = path_benign, family = binomial())

# Model summary
summary(logistic_model)

# Visualize results with a boxplot
boxplot(PHRED ~ pathogenic, data = path_benign, main = "PHRED by Pathogenicity")








### Weighted logistic regression 






prob_imputation <- function(data, prop_path){
  
  # Adding in weights 
  total_obs <- sum(!is.na(filtered_data$pathogenic)) #24176
  n_case <- sum(filtered_data$pathogenic == 1, na.rm=T) #439
  n_cont <- sum(filtered_data$pathogenic == -1, na.rm=T) #3364
  
  w_case <- total_obs * prop_path / n_case #2.202825
  w_cntl <- total_obs * (1 - prop_path) / n_cont #6.899215
  
  ## number of imputations used to come up with estimates of logistics regression 
  ## parameters.
  n_imputations = 25
  
  ## all data
  data_all <- data
  
  ## observed data (dis / benign) (remove VUS/uncertain which are listed as NAs)
  # Filter it so that pathogenic is either -1 (benign) or 1 (pathogenic) >> removing uncertain variants
  d_observed <- filtered_data %>% dplyr::filter(!pathogenic == 0) %>%
    dplyr::mutate(w = ifelse(pathogenic == 1, w_case, w_cntl)) #28481 observations
  
  ## uncertain data () containing the uncertain variants 
  d_unc <- filtered_data %>% dplyr::filter(pathogenic == 0) #should this be filtered to anything that = 0
  
  
  ## calculate AUC values and get data to make ROC curves
  #tmp <- AUCs <- auc_imputation(d_observed, d_unc)
  #AUCs <- tmp$aucs
  #roc_data <- tmp$roc_data
  
  model_estimates <- list()
  model_estimates_noweight <- list()
  model_estimates_weight <- list()
  plot_data <- c()
  
  #annot <- "PHRED"
  
  for (annot in c("PHRED")){
    # Doesn't need to be a loop since we're only using CADD (annot = CADD)
    #for (annot in c("CADD")){
    
    ## ############################
    ## Perform weighted logistic regression to get parameter estimate
    ## and estimate probability that uncertain variants are pathogenic
    ## ############################
    
    # Converting pathogenic variable to factor with levels -1 (benign) and 1 (pathogenic)
    d_obs <- d_observed %>% dplyr::select(all_of(c("merge", annot, "pathogenic", "w"))) %>%
      mutate(pathogenic=factor(pathogenic, levels = c("-1","1"))) 
    
    # Formula string where pathogenic is dependent variable and annot is independent variable
    formula <- paste("pathogenic ~", annot)
    
    # Converting the string to an actual formula
    formula <- as.formula(formula)
    
    # Adjust pathogenic variable to keep the 1 as 1 and change the -1 to 0
    d_obs <- d_obs %>%
      dplyr::mutate(pathogenic = as.numeric(ifelse(pathogenic == "-1", 0, as.character(pathogenic))))
    
    # Run logistic regression model with weights
    # Model is predicting "pathogenic" based on "annot" or PHRED
    model <- glm(formula, data = d_obs, family="binomial", weights=w) 
    #model <- glm(formula, data = d_obs, family="binomial")  
    
    
    # Can stay as binomial because there are no weights 
    # Run logistic regression without weights >> all observations are treated equally
    # Allows us to see how much the weights affect the model and use as means of comparison 
    model_noweight <- glm(formula, data = d_obs, family="binomial") 
    
    
    ## Estimate probably uncertain variants are pathogenic
    
    ## uncertain data () containing the uncertain variants 
    d_unc <- filtered_data %>% dplyr::filter(pathogenic == 0) #should this be filtered to anything that = 0
    
    # get rid of all columns besides merge, prob_disease)
    colnames(d_unc)
    d_unc <- d_unc %>% dplyr::select(merge, pathogenic, PHRED)
    
    # Use this to assign imputations a 1 or -1
    # Uses the logistic regression model to predict probability of being pathogenic for uncertain variants 
    # prob_disease is probability of pathogenic variant
    d_unc$prob_disease <- predict(model, newdata = d_unc, type = "response") #look at plot of pathogenic vs PHRED
    
    # Clean up d_unc again
    d_unc <- d_unc %>% dplyr::select(merge, prob_disease)
    
    
    ## reweight for the probabilistic imputation for uncertain variants that now have a probability
    # known benign and pathogenic get a weight of 1 (highest weight you can have bc we're confident)
    data_all <- filtered_data #replace with filtered_data 
    d_all <- data_all %>% dplyr::mutate(weight = 1) # creates new column called weight and sets all rows = 1
    d_all <- left_join(d_all, d_unc, by = "merge") #left join with filtered data
    # d_all <- d_all %>%
    #   mutate(
    #     pathogenic = coalesce(pathogenic.x, pathogenic.y),
    #     PHRED = coalesce(PHRED.x, PHRED.y)
    #   ) %>%
    #   dplyr::select(-pathogenic.x, -pathogenic.y, -PHRED.x, -PHRED.y)
    
    # Remove the old columns if needed (I didn't)
    
    
    #d_all <- left_join(d_all, d_unc, join_by(merge)) # dataframe with uncertain predictions -- will give prob_disease a value of NA if there isn't a value from the model
    #d_all <- d_all %>% dplyr::mutate(weight = ifelse(pathogenic.x == 0, abs(prob_disease - .5)*2, weight)) #change to pathogenic == 0
    
    # would be -0.5 or 0.5 if confident that variant is benign or pathogenic
    # want to see how close it is to 1 >> essentially turning it into a confidence
    d_all <- d_all %>%
      dplyr::mutate(weight = ifelse(pathogenic == 0, abs(prob_disease - 0.5) * 2, weight)) #0.5 is lowest weight you can have
    
    
    ## number of imputations used to come up with estimates of logistics regression parameters.
    n_imputations = 25
    models <- vector("list", n_imputations) # Store models
    aucs <- numeric(n_imputations) # Store AUC imputation
    
    
    ## perform naive estimates
    model_estimates_weight <- list(params = summary(model)$coef[,1],
                                   cov = vcov(model),
                                   t = summary(model)$coef[,3],
                                   p = summary(model)$coef[,4])
    # 
    # model_estimates_noweight[[annot]] <- list(params = summary(model_noweight)$coef[,1],
    #                                           cov = vcov(model_noweight), 
    #                                           t = summary(model_noweight)$coef[,3], 
    #                                           p = summary(model_noweight)$coef[,4])
    plot_data <- c()
    
    for (i in 1:n_imputations) {
      print(i)
      # Assuming 'df' has columns: 'predictor', 'classification' (NA for uncertain), 'prob_pathogenic'
      d_imp <- d_all
      
      # Assigning probabilistic classifications
      # Assigning pathogenic value where pathogenic = 0 
      # Going to randomly assign a 0 or 1 where pathogenic = 0
      # Sum of weights should be equal to number of pathogenic variants 
      # Results should be similar to original logistic regression 
      #### not finding pathogenic variable
      #d_imp <- d_imp %>% rowwise() %>%
      #dplyr::mutate(pathogenic = ifelse(pathogenic == 0, rbinom(1, 1, prob_disease), pathogenic)) #should we convert back to -1, 0, 1
      
      # reintroduce 0 and 1 (-1 turns to 0) 
      d_imp <- d_imp %>% 
        rowwise() %>%
        mutate(pathogenic = ifelse(!is.na(prob_disease), rbinom(1, 1, prob_disease), pathogenic),
               pathogenic = ifelse(pathogenic == -1, 0, pathogenic)) %>%
        ungroup()  # remove pathogenic == 0
      
      
      
      # Fit the weighted logistic regression model
      # Using above data to fit logistic regression with weights defined from confidence from pathogenic label
      models[[i]] <- glm(formula, family = binomial, data = d_imp, weights = weight)
      
      
      ## ###############
      ## Collect data for plotting
      # Storing data from first imputation
      ## ###############
      if (i == 1){
        tmp <- d_imp %>% dplyr::select(pathogenic, PHRED) 
        #dplyr::rename(score = !!sym(annot)) %>% 
        #dplyr::mutate(annot = annot)
        if (is.null(plot_data)){
          plot_data = tmp
        }else{
          plot_data = rbind(plot_data, tmp)
        }
      }
    }
    
    
    
    ## #####
    ## COMBINE INFORMATION ACCROSS IMPUTATIONS
    ## #####
    
    # Placeholder for aggregated variance-covariance matrices
    vcov_within <- matrix(0, nrow = length(models[[1]]$coefficients), ncol = length(models[[1]]$coefficients))
    vcov_between <- matrix(0, nrow = length(models[[1]]$coefficients), ncol = length(models[[1]]$coefficients))
    
    # Calculate W and B
    for(i in 1:n_imputations) {
      vcov_within <- vcov_within + vcov(models[[i]])
    }
    vcov_within <- vcov_within / n_imputations
    
    beta_bar <- sapply(models, function(x) coef(x)) %>% rowMeans
    for (i in 1:n_imputations) {
      diff <- matrix(coef(models[[i]]) - beta_bar, ncol = length(beta_bar))
      vcov_between <- vcov_between + (t(diff) %*% diff)
    }
    vcov_between <- vcov_between / (n_imputations - 1)
    
    # Calculate T
    vcov_total <- vcov_within + (1 + 1/length(models)) * vcov_between
    
    # T and P-value
    t_stat <- beta_bar / sqrt(diag(vcov_total))
    p_value <- 2 * pt(-abs(t_stat), df = length(models) - 1)
    
    model_estimates[[annot]] <- list(params = beta_bar, cov = vcov_total, t = t_stat, p = p_value)
    
    
  }
  
  
  # model_estimates will then be used to create the table
  
  return(list(model_estimates = model_estimates, 
              #model_estimates_noweight = model_estimates_noweight, 
              model_estimates_weight = model_estimates_weight, 
              aucs = AUCs, 
              roc_data = roc_data, 
              plot_data = plot_data))
}



## #################
## Sneha Threshold Table 
## ################

### STEP 1: Define Pathogenicity Proportion for Weighting
prop_path_update <- 0.04  # Estimated proportion of pathogenic variants

### STEP 2: Run Weighted Logistic Regression with Imputation
tmp <- prob_imputation(filtered_data, prop_path = prop_path_update)

### STEP 3: Extract Model Results
est_perm <- tmp$model_estimates  # Required for ACMG threshold calculation

est_noweight <- tmp$model_estimates_noweight  # Optional (comparison without weights)
est_weight <- tmp$model_estimates_weight  # Optional (comparison with weights)
aucs <- tmp$aucs  # Optional (check model performance)
roc_data <- tmp$roc_data  # Optional (ROC curve data)
plot_data <- tmp$plot_data  # Optional (Pathogenicity visualization)

### STEP 4: Compute ACMG Thresholds
calculate_thresholds <- function(est_perm, filtered_data, prop_path, file="acmg_table_impute.txt") {
  
  rslts <- c()
  cutoffs <- c(2.406, 5.790, 33.53, 1124)  # ACMG cutoffs
  cutlabel <- c("supporting", "moderate", "strong", "very_strong") 
  pD <- prop_path
  odds_path <- pD / (1 - pD)
  
  # Function to compute pathogenicity thresholds
  min_func <- function(x, beta, v, cutoff, odds_path, var_type) {
    if (!var_type %in% c("path", "benign")) {
      stop("var_type must be 'path' or 'benign'")
    }
    
    var_adj <- ifelse(var_type == "path", 1, -1)  # Adjust based on pathogenic or benign
    se_log_odds <- sqrt(t(c(1, x)) %*% v %*% c(1, x))  # Compute standard error
    z <- qnorm(.95)  # 95% CI
    logodds_bound <- (beta[1] + beta[2] * x) - var_adj * z * se_log_odds
    odds_bound <- exp(var_adj * logodds_bound)
    odds_path <- odds_path^var_adj
    
    return(cutoff - odds_bound / odds_path)
  }
  
  # Iterate over annotation scores (only PHRED in your case)
  for (annot in c("PHRED")) {
    
    betas <- est_perm[[annot]]$params
    v <- est_perm[[annot]]$cov
    int <- betas[1]   # Intercept
    logor <- betas[2]  # Log odds ratio
    p <- est_perm[[annot]]$p[2] # P-value
    se_int <- sqrt(v[1, 1])
    se_logor <- sqrt(v[2, 2])
    
    rng <- range(filtered_data[[annot]], na.rm = TRUE)  # Find PHRED range
    lower <- min(rng) 
    upper <- max(rng)  
    
    # Loop through ACMG categories
    for (j in seq_along(cutoffs)) {
      cutoff <- cutoffs[j]
      
      # Compute pathogenic threshold
      result <- tryCatch({
        uniroot(min_func, lower = lower, upper = upper, beta = betas, v = v,  
                cutoff = cutoff, odds_path = odds_path, var_type = "path")
      }, error = function(e) list(root = NA))
      
      rslts <- rbind(rslts, c(annot, int, logor, se_int, se_logor, p, 
                              "pathogenic", cutlabel[j], result$root, rng))
      
      # Compute benign threshold
      result <- tryCatch({
        uniroot(min_func, lower = lower, upper = upper, beta = betas, v = v,
                cutoff = cutoff, odds_path = odds_path, var_type = "benign")
      }, error = function(e) list(root = NA))
      
      rslts <- rbind(rslts, c(annot, int, logor, se_int, se_logor, p, 
                              "benign", cutlabel[j], result$root, rng))
    }       
  }
  
  # Convert to dataframe
  rslts <- as.data.frame(rslts, stringsAsFactors = FALSE)
  names(rslts) <- c("annot", "int", "logor", "se_int", "se_logor", "pval", "var_type", 
                    "acmg_cat", "annot_value", "min_obs", "max_obs")
  
  # Reshape for readability
  tbl <- rslts %>% tidyr::pivot_wider(names_from = c("var_type", "acmg_cat"),
                                      values_from = "annot_value",
                                      names_glue = "{var_type}_{acmg_cat}") %>%
    dplyr::select(annot, int, logor, se_int, se_logor, pval, 
                  starts_with("pathogenic"), starts_with("benign"),
                  min_obs, max_obs)
  
  # Save table to file
  write.table(file = file, tbl, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  cat("Wrote ACMG threshold table to", file, "\n")
  
  return(tbl)
}

### STEP 5: Run ACMG Threshold Calculation
threshold_table <- calculate_thresholds(est_perm, filtered_data, prop_path_update)

print(threshold_table)







### STOPPED HERE -- EVERYTHING ABOVE HERE WORKS AS EXPECTED





## #########
## Write auc table ##
## #########
write.table(file=paste0(results_dir, "AUC_imputed_80_20.txt"), 
            roc_data, quote=F, row.names=F, col.names=T, sep="\t")

## ##############
## plot ROC curve using an imputed dataset ##
## ##############
plot_file <- paste0(results_dir,"FullSample_ROC_plot.pdf")
plot_roc_imputed(plot_data, plotfile = plot_file)

tbl <- acmg_table_impute(est_perm, data=data, prop_path = prop_path_update,
                         file=paste0(results_dir,"acmg_table_impute.txt"))
















auc_imputation <- function(d_obs, d_unc, annot, p_train=.8, n_imp=200){
  
  AUCs <- c()
  roc_data = c()
  aucs <- numeric(n_imp) # Store AUC imputation
  n_path = sum(d_obs$dis_ben == 1, na.rm=T)
  n_ben = sum(d_obs$dis_ben == 0, na.rm=T)
  n_unc = nrow(d_unc)
  
  for (annot in c("CADD","REVEL","BayesDel", "AlphaMissense")){
    
    cat("Working on ", annot,"\n")
    
    for (i in 1:n_imp){
      
      i_path <- sample(1:n_path,  ceiling(n_path * p_train))
      i_ben  <- sample(1:n_ben,   ceiling(n_ben * p_train))
      i_unc  <- sample(1:n_unc,   ceiling(n_unc * p_train))
      
      ## ##########
      ## Create sampled tested and training datasets
      ## ##########
      
      ## training sets ##
      d_path_tr = d_obs %>% dplyr::filter(dis_ben == 1) %>% 
        dplyr::mutate(row = row_number()) %>%
        dplyr::filter(row %in% i_path)
      d_ben_tr = d_obs %>% dplyr::filter(dis_ben == 0) %>% 
        dplyr::mutate(row = row_number()) %>%
        dplyr::filter(row %in% i_ben)
      d_unc_tr <- d_unc %>% dplyr::mutate(row = row_number()) %>%
        dplyr::filter(row %in% i_unc)
      
      ## testing sets
      d_path_te = d_obs %>% dplyr::filter(dis_ben == 1) %>% 
        dplyr::mutate(row = row_number()) %>%
        dplyr::filter(row %in% i_path == F)
      d_ben_te = d_obs %>% dplyr::filter(dis_ben == 0) %>% 
        dplyr::mutate(row = row_number()) %>%
        dplyr::filter(row %in% i_ben == F)
      d_unc_te <- d_unc %>% dplyr::mutate(row = row_number()) %>%
        dplyr::filter(row %in% i_unc == F)
      
      ## ###############
      ## Add the prob of disease for uncertain variants
      ## ################
      
      ### Get weights for unc variants ##
      ## Training
      
      formula <- paste("dis_ben ~", annot)
      # Converting the string to an actual formula
      formula <- as.formula(formula)
      
      d_o_tr <- rbind(d_path_tr, d_ben_tr) %>%
        dplyr::select(all_of(c(annot, "dis_ben", "w")))
      model <- glm(formula, data = d_o_tr, family="binomial", weight=w)
      
      d_unc_tr$prob_disease <- predict(model, newdata = d_unc_tr, type = "response") 
      d_unc_tr <- d_unc_tr %>% 
        rowwise() %>%
        dplyr::mutate(dis_ben = rbinom(1, 1, prob_disease), 
                      w = abs(prob_disease - .5)*2 ) %>%
        dplyr::select(all_of(c(annot, "dis_ben", "w"))) %>% ungroup()               
      
      d_o_te <- rbind(d_path_te, d_ben_te) %>%
        dplyr::select(all_of(c(annot, "dis_ben", "w")))
      model <- glm(formula, data = d_o_te, family="binomial", weight=w)
      d_unc_te$prob_disease <- predict(model, newdata = d_unc_te, type = "response")
      d_unc_te <- d_unc_te %>% 
        rowwise() %>%
        dplyr::mutate(dis_ben = rbinom(1, 1, prob_disease), 
                      w = abs(prob_disease - .5)*2 ) %>%
        dplyr::select(all_of(c(annot, "dis_ben", "w")))       
      
      ## ###########
      ## Reweight observed data for the probabilistic imputation 
      ## Combine observed and uncertain
      ## #########
      d_o_tr$w = 1; d_o_te$w = 1
      d_all_tr <- rbind(d_o_tr, d_unc_tr)
      d_all_te <- rbind(d_o_te, d_unc_te)
      
      ## ########
      ## Fit the weighted logistic regression model to test sample
      ## ########
      model <- glm(formula, family = binomial, data = d_all_tr, weights = w)
      preds <- predict(model, d_all_te, type="response")
      
      ## ####
      ## Calculate AUC 
      ## ####
      pred_obj <- ROCR::prediction(preds, d_all_te$dis_ben)  # actual_responses needs to match the model's data
      
      auc <- performance(pred_obj, "auc")
      aucs[i] <- auc@y.values[[1]]  # Extract AUC value
      
      ## ###############
      ## Collect data for plotting 
      ## ###############
      if (i == 1){
        tmp <- d_all_te %>% dplyr::select(all_of(c("dis_ben", annot))) %>%
          dplyr::rename(score = !!sym(annot)) %>% 
          dplyr::mutate(annot = annot)
        if (is.null(roc_data)){
          roc_data = tmp
        }else{
          roc_data = rbind(roc_data, tmp)
        }
      }
    }
    AUCs <- rbind(AUCs, c(annot, mean(aucs), sqrt(var(aucs))))
  }
  
  return(list(aucs = AUCs, roc_data = roc_data))
}

#Threshold Table

## #################
<<<<<<< HEAD
## INCORPORATE UNCERTAIN VARIANTS USING IMPUTATION
## ################
prop_path_update <- gnomad$prob_path_update # estimate of proportion of pathogenic variants (just use 0.04 for now) 
prop_path_update <- 0.04
tmp <- prob_imputation(data, prop_path = prop_path_update)

est_perm <- tmp$model_estimates
est_noweight <- tmp$model_estimates_noweight
est_weight <- tmp$model_estimates_weight
aucs <- tmp$aucs
roc_data <- tmp$roc_data
plot_data <- tmp$plot_data

acmg_table_impute <- function(est_perm, data, prop_path, file="acmg_table_impute.txt"){
  
  rslts <- c()
  cutoffs <- c(2.406, 5.790, 33.53, 1124)
  cutlabel <- c("supporting", "moderate", "strong", "very_strong")
  pD = prop_path
  odds_path = pD / (1 - pD)
  
  
  min_func <- function(x, beta, v, cutoff, odds_path, var_type ){
    
    if (! var_type %in% c("path","benign")){
      cat("var_type must be path or benign\n")
      stop()
    }
    
    ## calc se of log odds ##
    var_adj <- ifelse(var_type == "path", 1, -1)
    se_log_odds <- sqrt(t(c(1, x)) %*% v %*% c(1, x)) # Standard error
    z <- qnorm(.95)
    logodds_bound <- (beta[1] + beta[2]*x) - var_adj * z * se_log_odds
    odds_bound <- exp( var_adj * logodds_bound)
    odds_path = odds_path^var_adj
=======
## Sneha Threshold Table 
## ################
    
### STEP 1: Define Pathogenicity Proportion for Weighting
prop_path_update <- 0.04  # Estimated proportion of pathogenic variants

### STEP 2: Run Weighted Logistic Regression with Imputation
tmp <- prob_imputation(filtered_data, prop_path = prop_path_update)

### STEP 3: Extract Model Results
est_perm <- tmp$model_estimates  # Required for ACMG threshold calculation

est_noweight <- tmp$model_estimates_noweight  # Optional (comparison without weights)
est_weight <- tmp$model_estimates_weight  # Optional (comparison with weights)
aucs <- tmp$aucs  # Optional (check model performance)
roc_data <- tmp$roc_data  # Optional (ROC curve data)
plot_data <- tmp$plot_data  # Optional (Pathogenicity visualization)

### STEP 4: Compute ACMG Thresholds
calculate_thresholds <- function(est_perm, filtered_data, prop_path, file="acmg_table_impute.txt") {
  
  rslts <- c()
  cutoffs <- c(2.406, 5.790, 33.53, 1124)  # ACMG cutoffs
  cutlabel <- c("supporting", "moderate", "strong", "very_strong") 
  pD <- prop_path
  odds_path <- pD / (1 - pD)
  
  # Function to compute pathogenicity thresholds
  min_func <- function(x, beta, v, cutoff, odds_path, var_type) {
    if (!var_type %in% c("path", "benign")) {
      stop("var_type must be 'path' or 'benign'")
    }
    
    var_adj <- ifelse(var_type == "path", 1, -1)  # Adjust based on pathogenic or benign
    se_log_odds <- sqrt(t(c(1, x)) %*% v %*% c(1, x))  # Compute standard error
    z <- qnorm(.95)  # 95% CI
    logodds_bound <- (beta[1] + beta[2] * x) - var_adj * z * se_log_odds
    odds_bound <- exp(var_adj * logodds_bound)
    odds_path <- odds_path^var_adj
>>>>>>> fbc16e7b564b02a8abbf47a818e2f3b5b9e65019
    
    return(cutoff - odds_bound / odds_path)
  }
  
<<<<<<< HEAD
  for (i in 1:4){
    
    betas = est_perm[[i]]$params
    v = est_perm[[i]]$cov
    cat = names(est_perm)[i]
    int <- betas[1]
    logor <- betas[2] 
    p <- est_perm[[i]]$p[2]
    se_int <- sqrt(v[1,1])
    se_logor <- sqrt(v[2,2])
    
    rng <- range(data[, cat])
    lower = -50; upper = 50
    #if (i == 4){ lower = 0; upper=1} ## alphamissense
    
    for (j in 1:length(cutoffs)){
      cutoff <- cutoffs[j]
      
      #minimize_func_path <- function(x){ return( cutoff - exp(int + x*lor)/odds_path) }
      #minimize_func_ben <- function(x){ return( cutoff -  exp(-(int + x*lor))*odds_path) }
      
      result <- tryCatch({
        uniroot(min_func, lower = lower, upper = upper, beta=betas, v = v,  
                cutoff = cutoff, odds_path = odds_path, var_type = "path")
      }, error = function(e){
        cat("No satisfying threshold was able to be found for:\n")
        cat(cat, " using cutoff: ", cutoff, " for pathogenic variant\n")
        return(list(root = NA))
      })
      
      rslts <- rbind(rslts, c(cat, int, logor, se_int, se_logor, p, 
                              "pathogenic", cutlabel[j], result$root, rng))
      
      result <- tryCatch({
        uniroot(min_func, lower = lower, upper = upper, beta=betas, v = v,
                cutoff = cutoff, odds_path = odds_path, var_type = "benign")
      }, error = function(e){
        cat("No satisfying threshold was able to be found for:\n")
        cat(cat, " using cutoff: ", cutoff, " for benign variant\n")
        return(list(root = NA))
      })
      rslts <- rbind(rslts, c(cat, int, logor, se_int, se_logor, p, 
=======
  # Iterate over annotation scores (only PHRED for now)
  for (annot in c("PHRED")) {
    
    betas <- est_perm[[annot]]$params
    v <- est_perm[[annot]]$cov
    int <- betas[1]   # Intercept
    logor <- betas[2]  # Log odds ratio
    p <- est_perm[[annot]]$p[2] # P-value
    se_int <- sqrt(v[1, 1])
    se_logor <- sqrt(v[2, 2])
    
    rng <- range(filtered_data[[annot]], na.rm = TRUE)  # Find PHRED range
    lower <- min(rng) 
    upper <- max(rng)  
    
    # Loop through ACMG categories
    for (j in seq_along(cutoffs)) {
      cutoff <- cutoffs[j]
      
      # Compute pathogenic threshold
      result <- tryCatch({
        uniroot(min_func, lower = lower, upper = upper, beta = betas, v = v,  
                cutoff = cutoff, odds_path = odds_path, var_type = "path")
      }, error = function(e) list(root = NA))
      
      rslts <- rbind(rslts, c(annot, int, logor, se_int, se_logor, p, 
                              "pathogenic", cutlabel[j], result$root, rng))
      
      # Compute benign threshold
      result <- tryCatch({
        uniroot(min_func, lower = lower, upper = upper, beta = betas, v = v,
                cutoff = cutoff, odds_path = odds_path, var_type = "benign")
      }, error = function(e) list(root = NA))
      
      rslts <- rbind(rslts, c(annot, int, logor, se_int, se_logor, p, 
>>>>>>> fbc16e7b564b02a8abbf47a818e2f3b5b9e65019
                              "benign", cutlabel[j], result$root, rng))
    }       
  }
  
<<<<<<< HEAD
  rslts <- as.data.frame(rslts, stringsAsFactors=F)
  names(rslts) <- c("annot", "int", "logor", "se_int", "se_logor", "pval", "var_type", 
                    "acmg_cat", "annot_value", "min_obs", "max_obs")
  
  tbl <- rslts %>% tidyr::pivot_wider(names_from=c("var_type","acmg_cat"),
                                      values_from="annot_value",
=======
  # Convert to dataframe
  rslts <- as.data.frame(rslts, stringsAsFactors = FALSE)
  names(rslts) <- c("annot", "int", "logor", "se_int", "se_logor", "pval", "var_type", 
                    "acmg_cat", "annot_value", "min_obs", "max_obs")
  
  # Reshape for readability
  tbl <- rslts %>% tidyr::pivot_wider(names_from = c("var_type", "acmg_cat"),
                                      values_from = "annot_value",
>>>>>>> fbc16e7b564b02a8abbf47a818e2f3b5b9e65019
                                      names_glue = "{var_type}_{acmg_cat}") %>%
    dplyr::select(annot, int, logor, se_int, se_logor, pval, 
                  starts_with("pathogenic"), starts_with("benign"),
                  min_obs, max_obs)
  
<<<<<<< HEAD
  
  write.table(file = file, tbl, quote=F, row.names=F, col.names=T, sep="\t")
  cat("Wrote threshold for different ACMG support categories to", file, "\n")
=======
  # Save table to file
  write.table(file = file, tbl, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  cat("Wrote ACMG threshold table to", file, "\n")
>>>>>>> fbc16e7b564b02a8abbf47a818e2f3b5b9e65019
  
  return(tbl)
}

<<<<<<< HEAD
=======
### STEP 5: Run ACMG Threshold Calculation
threshold_table <- calculate_thresholds(est_perm, filtered_data, prop_path_update)

print(threshold_table)
>>>>>>> fbc16e7b564b02a8abbf47a818e2f3b5b9e65019

## #################
## CI of Threshold Table 
## ################

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Define number of imputations
n_imputations <- 25  # Number of imputed datasets

# Store model results across imputations
beta_values <- numeric(n_imputations)
se_values <- numeric(n_imputations)

# Ensure prob_disease is calculated from Model 1
filtered_data <- filtered_data %>%
  mutate(prob_disease = predict(logistic_model, newdata = filtered_data, type = "response"))

# Ensure weight column is correctly assigned
filtered_data <- filtered_data %>%
  mutate(weight = ifelse(pathogenic == 1, 1, ifelse(pathogenic == -1, 1, abs(prob_disease - 0.5) * 2)))

# Loop over imputations, running logistic regression for each one
for (i in 1:n_imputations) {
  
  # Create a new imputed dataset by assigning probabilistic classifications to VUS
  d_imp <- filtered_data %>%
    rowwise() %>%
    mutate(pathogenic = ifelse(!is.na(prob_disease), rbinom(1, 1, prob_disease), pathogenic),
           pathogenic = ifelse(pathogenic == -1, 0, pathogenic)) %>%
    ungroup()
  
  # Ensure weight column is present in imputed dataset
  d_imp <- d_imp %>%
    mutate(weight = ifelse(pathogenic == 1, 1, ifelse(pathogenic == 0, abs(prob_disease - 0.5) * 2, 1)))
  
  # Fit weighted logistic regression
  model <- glm(pathogenic ~ PHRED, family = binomial, data = d_imp, weights = weight)
  
  # Extract beta coefficient and standard error for PHRED
  beta_values[i] <- coef(model)["PHRED"]
  se_values[i] <- summary(model)$coefficients["PHRED", "Std. Error"]
}

# Calculate the mean beta coefficient across imputations
beta_avg <- mean(beta_values)

# Calculate within-imputation variance (W)
var_within <- mean(se_values^2)

# Calculate between-imputation variance (B)
var_between <- var(beta_values)

# Total variance using Rubin’s rules
var_total <- var_within + (1 + (1/n_imputations)) * var_between

# Compute 95% Confidence Interval
z_value <- 1.96  # For 95% CI
lower_CI <- beta_avg - z_value * sqrt(var_total)
upper_CI <- beta_avg + z_value * sqrt(var_total)

# Print results
cat("Beta Estimate for CADD PHRED (Weighted Model):", beta_avg, "\n")
cat("95% Confidence Interval: (", lower_CI, ",", upper_CI, ")\n")

# Plot the beta estimates across imputations
ggplot(data.frame(beta_values), aes(x = beta_values)) +
  geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.6, color = "black") +
  labs(title = "Distribution of Beta Estimates Across Imputations",
       x = "Beta Estimate", y = "Frequency") +
  theme_minimal()




## #################
## Table/ Results Ideas
## ################

# Summary Statistics: Table 1

summary_table <- data.frame(
  Source = c("ClinVar", "gnomAD Benign"),
  Variants_Before_Filtering = c(nrow(pten_df), nrow(gnomad_fixed)),  # Raw counts before merging
  Variants_After_Filtering = c(
    nrow(final_combined_data_2[!is.na(final_combined_data_2$CLNSIG),]),  # ClinVar variants after filtering
    nrow(final_combined_data_2[is.na(final_combined_data_2$CLNSIG),])   # gnomAD variants after filtering
  ),
  Pathogenic = c(sum(final_combined_data_2$CLNSIG %in% c("Pathogenic", "Likely_pathogenic"), na.rm=TRUE), NA),
  Benign = c(sum(final_combined_data_2$CLNSIG %in% c("Benign", "Likely_benign"), na.rm=TRUE), NA),
  Uncertain = c(sum(final_combined_data_2$CLNSIG == "Uncertain", na.rm=TRUE), NA)
)

print(summary_table)





## #################
## Ideas to show results
## ################

# Heat Map for Thresholds 

install.packages("reshape2")
library(reshape2)
library(ggplot2)

# Convert threshold table values to numeric (except "annot" column)
threshold_melt <- melt(threshold_table, id.vars = "annot")
threshold_melt$value <- as.numeric(as.character(threshold_melt$value))

# Plot heatmap
ggplot(threshold_melt, aes(x = variable, y = annot, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "ACMG Pathogenicity Thresholds", x = "ACMG Category", y = "Annotation Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Logistic Regression Model Preformance
summary(tmp$model_estimates$PHRED)

roc_curve <- roc(path_benign$pathogenic, predict(tmp$model_estimates$PHRED[[1]], newdata = path_benign, type="response"))

ggplot(data = data.frame(FPR = roc_curve$specificities, TPR = roc_curve$sensitivities), 
       aes(x = FPR, y = TPR)) +
  geom_line(color = "blue") +
  geom_abline(linetype="dashed", color="red") +
  labs(title = "ROC Curve for Pathogenicity Prediction",
       x = "False Positive Rate", y = "True Positive Rate")

auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))






### Annotation score as a function of position plot

# Exon positions
exonfile <- '~/GitHub/Capstone_2024/Resources/CodingLocs_hg38.rds'

readRDS(exonfile)
exonpositions <- readRDS(exonfile)

















# Plot
homeo_plot <- function(data_all, prefix = "homeo"){
  
  ## plot annotation scores as a function of AA position
  
  plot_func <- function(d){
    
    ## values from Andy via Alamute Visual.
    hd <- c(292, 471) #amino acid position (homeodomain)
    exons <- c(241, 429, 945) #UCSC file with genomic positions of start and stop of every axon
      # adjust the positions so that first exon is end - start and second one is first positoin +1 = total number of base pairs in eawhc exon (write a function for it)
      # also plot boundaries of exons as vertical lines
      # would be worth seeing if variants fall into exos
     p <- ggplot(d, aes(x = pos_c, y = score, color = pathannot))
    
    p <- p + geom_point(size=2, alpha=.4) +
      facet_wrap(~factor(annot), scales = "free_y") 
    p <- p + scale_color_manual(values = c("Benign" = "darkgreen", "Pathogenic"="darkred", 
                                           "Uncertain" = "blue"))
    p <- p + geom_vline(xintercept = hd, linetype = "dashed", color = "red")
    p <- p + geom_vline(xintercept = exons)
    p <- p + theme_bw() + xlab("Position") 
    p <- p + labs(color = "Variant Class") # Change title for color legend)
    p <- p + guides(color = guide_legend(override.aes = list(shape = 16 ))) + # Remove characters from color legend
      theme(
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), # Hide x-axis ticks and text
        strip.text = element_text(size = 14, face = "bold"), # Facet label size and bold
        axis.title.y = element_text(size = 14, face = "bold", ) # Y-axis title size and bold
      )
    
    file <- paste0(prefix, "_score_x_pos.pdf")
    
    pdf(file, width=11, height = 8.5, )
    print(p)
    dev.off()
    cat("Created plot", file, "\n")
  }
  
  scores <- colnames(data_all)[colnames(data_all) %in% c("PHRED")]
  
  d <- data_all %>% dplyr::mutate(exon = factor(exon, levels = c("1", "2", "3")))
  d <- d %>% tidyr::pivot_longer(cols = all_of(scores), names_to = "annot", values_to = "score") %>%
    dplyr::mutate(annot = factor(annot, levels = c("PHRED")))
  
  plot_func(d)
}

# Exon function to integrate above (figure out what the positions we're referring to [p., c., etc])
# https://genome.ucsc.edu/cgi-bin/hgTables
# In essence, this function translates a genomic coordinate into a coding sequence coordinate, taking into account the exon structure of a gene.
gene_info <- read.table(file = "~/GitHub/ClinVar-GnomAD-Merge/pten_exon_postions", as.is=T, header=F, sep="\t")
names(gene_info)[1:3] <- c("chr", "start", "end")

find_cdot <- function(pos, gene_info){
  exon_size <- gene_info$end - gene_info$start
  pos_adj <- cumsum(exon_size)
  pos_adj <- c(0, pos_adj[-length(pos_adj)]) #shifts the cumulative sums one position to the right and adds a 0 at the beginning. This creates a vector where each element represents the starting position of the next exon in the coding sequence.
  gene_info$exon <- 1:nrow(gene_info) #This line adds an exon column to the gene_info table, assigning sequential numbers from 1 to the number of rows (exons).
  result <- gene_info %>% 
    dplyr::mutate(in_here = start <= pos  & end >= pos,
                  pos_wi_exon = pos - start) %>% 
    dplyr::filter(in_here == TRUE) %>%
    dplyr::select(exon, pos_wi_exon) 
  coding_pos <- pos_adj[result$exon] + result$pos_wi_exon + 1 
  return(coding_pos)
}


test_pos <- 87958019 - 30
find_cdot(test_pos, gene_info)













# Load libraries (if not already loaded)
library(dplyr)

# Assuming d_all and gene_info are already loaded

# Add exon column to gene_info (if not already present)
if (!("exon" %in% colnames(gene_info))) {
  gene_info$exon <- 1:nrow(gene_info)
}

# find_cdot function (as previously defined)
find_cdot <- function(pos, gene_info) {
  exon_size <- gene_info$end - gene_info$start
  pos_adj <- cumsum(exon_size)
  pos_adj <- c(0, pos_adj[-length(pos_adj)])
  gene_info$exon <- 1:nrow(gene_info)
  result <- gene_info %>%
    dplyr::mutate(in_here = start <= pos & end >= pos,
                  pos_wi_exon = pos - start) %>%
    dplyr::filter(in_here == TRUE) %>%
    dplyr::select(exon, pos_wi_exon)
  
  if (nrow(result) > 0) {
    coding_pos <- pos_adj[result$exon] + result$pos_wi_exon + 1
    return(coding_pos)
  } else {
    return(NA) # Return NA when no exon is found
  }
}

# Apply find_cdot to d_all$start
d_all$coding_start <- sapply(d_all$start.x, function(pos) find_cdot(pos, gene_info))

# Apply find_cdot to d_all$end
d_all$coding_end <- sapply(d_all$end.x, function(pos) find_cdot(pos, gene_info))

# Print the results
print(d_all)



# Add exon column to d_all (if not already present)
if (!("exon" %in% colnames(d_all))) {
  d_all$exon <- sapply(d_all$start, function(pos) {
    gene_info$exon[gene_info$start <= pos & gene_info$end >= pos]
  })
}

print(head(d_all)) # Verify exon before homeo_plot

homeo_plot <- function(data_all, prefix = "homeo", gene_info, remove_na = TRUE) {
  
  plot_func <- function(d) {
    exons <- gene_info$start[gene_info$exon != 1 & !is.na(gene_info$start)]
    p <- ggplot(d, aes(x = start.x, y = score, color = pathogenic)) #using start.x
    p <- p + geom_point(size = 2, alpha = 0.4) + facet_wrap(~factor(annot), scales = "free_y")
    p <- p + scale_color_manual(values = c("Benign" = "darkgreen", "Pathogenic" = "darkred", "Uncertain" = "blue"))
    if(length(exons) > 0){ p <- p + geom_vline(xintercept = exons, linetype = "dotted", color = "darkgray") }
    p <- p + theme_bw() + xlab("Position") + labs(color = "Variant Class")
    p <- p + guides(color = guide_legend(override.aes = list(shape = 16))) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 14, face = "bold"))
    file <- paste0(prefix, "_score_x_pos.pdf"); pdf(file, width = 11, height = 8.5); print(p); dev.off(); cat("Created plot", file, "\n")
  }
  
  scores <- colnames(data_all)[colnames(data_all) %in% c("PHRED")]
  
  print(table(is.na(data_all$coding_start)))
  print(table(is.na(data_all$exon)))
  print(paste("d_all row count before coding_start filter:", nrow(data_all)))
  
  if(remove_na == TRUE){
    data_all <- data_all %>% dplyr::filter(!is.na(coding_start))
  }
  
  print(paste("d_all row count after coding_start filter:", nrow(data_all)))
  
  d <- data_all %>%
    dplyr::filter(!is.na(exon)) %>%
    dplyr::mutate(exon = factor(exon, levels = as.character(unique(gene_info$exon))),
                  pathogenic = factor(case_when(pathogenic == 0 ~ "Benign", pathogenic == 1 ~ "Pathogenic", TRUE ~ "Uncertain"), levels = c("Benign", "Pathogenic", "Uncertain"))
    ) %>%
    tidyr::pivot_longer(cols = all_of(scores), names_to = "annot", values_to = "score") %>%
    dplyr::mutate(annot = factor(annot, levels = c("PHRED")))
  
  print(paste("d row count:", nrow(d)))
  print(head(d))
  print(table(d$start.x))
  print(sum(duplicated(d)))
  
  print(paste("d_all row count:", nrow(data_all)))
  print(paste("d row count:", nrow(d)))
  
  print(head(data.frame(data_all$coding_start, d$start.x)))
  
  plot_func(d)
}

if (!("exon" %in% colnames(d_all))) { d_all$exon <- sapply(d_all$start, function(pos) { gene_info$exon[gene_info$start <= pos & gene_info$end >= pos] }) }
homeo_plot(data_all = d_all, prefix = "my_d_all_plot", gene_info = gene_info, remove_na = FALSE)







# Calculating c. according to Gemini and plotting it 
calculate_coding_position <- function(variant_start, gene_info) {
  # Find the exon the variant is in
  exon_row <- gene_info[gene_info$start <= variant_start & gene_info$end >= variant_start, ]
  
  # If not in an exon, return NA
  if (nrow(exon_row) == 0) {
    return(NA)
  }
  
  exon_number <- exon_row$exon
  
  # Calculate offset from exon start
  offset <- variant_start - exon_row$start + 1
  
  # Calculate coding position
  coding_position <- 0
  
  # Sum lengths of upstream exons
  upstream_exons <- gene_info[gene_info$exon < exon_number, ]
  if (nrow(upstream_exons) > 0) {
    coding_position <- sum(upstream_exons$end - upstream_exons$start + 1)
  }
  
  # Add offset
  coding_position <- coding_position + offset
  
  return(coding_position)
}

homeo_plot <- function(data_all, prefix = "homeo", gene_info, remove_na = TRUE) {
  
  data_all$coding_start <- sapply(data_all$start.x, function(x) calculate_coding_position(x, gene_info))
  
  print(head(data_all$coding_start)) #checking coding_start
  
  plot_func <- function(d) {
    exons <- gene_info$start[gene_info$exon != 1 & !is.na(gene_info$start)]
    p <- ggplot(d, aes(x = pos_c, y = score, color = pathogenic))
    p <- p + geom_point(size = 2, alpha = 0.4) + facet_wrap(~factor(annot), scales = "free_y")
    p <- p + scale_color_manual(values = c("Benign" = "darkgreen", "Pathogenic" = "darkred", "Uncertain" = "blue"))
    if(length(exons) > 0){ p <- p + geom_vline(xintercept = exons, linetype = "dotted", color = "darkgray") }
    p <- p + theme_bw() + xlab("Position") + labs(color = "Variant Class")
    p <- p + guides(color = guide_legend(override.aes = list(shape = 16))) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 14, face = "bold"))
    file <- paste0(prefix, "_score_x_pos.pdf"); pdf(file, width = 11, height = 8.5); print(p); dev.off(); cat("Created plot", file, "\n")
  }
  
  scores <- colnames(data_all)[colnames(data_all) %in% c("PHRED")]
  
  print(table(is.na(data_all$coding_start)))
  print(table(is.na(data_all$exon)))
  print(paste("d_all row count before coding_start filter:", nrow(data_all)))
  
  # Filter NA in coding_start before creating 'd'
  data_all <- data_all %>% dplyr::filter(!is.na(coding_start))
  
  print(paste("d_all row count after coding_start filter:", nrow(data_all)))
  
  d <- data_all %>%
    dplyr::filter(!is.na(exon)) %>%
    dplyr::rename(pos_c = coding_start) %>%
    dplyr::mutate(exon = factor(exon, levels = as.character(unique(gene_info$exon))),
                  pathogenic = factor(case_when(pathogenic == 0 ~ "Benign", pathogenic == 1 ~ "Pathogenic", TRUE ~ "Uncertain"), levels = c("Benign", "Pathogenic", "Uncertain"))
    ) %>%
    tidyr::pivot_longer(cols = all_of(scores), names_to = "annot", values_to = "score") %>%
    dplyr::mutate(annot = factor(annot, levels = c("PHRED")))
  
  print(paste("d row count:", nrow(d)))
  print(head(d))
  print(table(d$pos_c))
  print(sum(duplicated(d)))
  
  print(paste("d_all row count:", nrow(data_all)))
  print(paste("d row count:", nrow(d)))
  
  print(head(data.frame(data_all$coding_start, d$pos_c)))
  
  print(gene_info) #checking gene_info
  print(head(data_all$start.x)) #checking start.x
  
  plot_func(d)
}

if (!("exon" %in% colnames(d_all))) { d_all$exon <- sapply(d_all$start, function(pos) { gene_info$exon[gene_info$start <= pos & gene_info$end >= pos] }) }
homeo_plot(data_all = d_all, prefix = "my_d_all_plot", gene_info = gene_info, remove_na = FALSE)



# Plotting according to Chatgpt -- this is the one that works 
library(dplyr)
library(ggplot2)
library(purrr)  # Load purrr for map_dbl() and map_chr()


# Load gene info
gene_info <- read.table(file = "~/GitHub/ClinVar-GnomAD-Merge/pten_exon_postions", 
                        as.is = TRUE, header = FALSE, sep = "\t")
names(gene_info)[1:3] <- c("chr", "start", "end")
gene_info$exon <- 1:nrow(gene_info)

# Compute cumulative exon lengths
gene_info <- gene_info %>%
  mutate(exon_length = end - start + 1,
         cumulative_length = cumsum(exon_length) - exon_length)

# Function to compute c. position (using start and end)
find_cdot <- function(start_pos, end_pos) {
  pos <- ifelse(is.na(end_pos) | start_pos == end_pos, start_pos, floor((start_pos + end_pos) / 2))
  exon_found <- gene_info %>%
    filter(start <= pos & end >= pos) %>%
    mutate(pos_within_exon = pos - start + 1)
  
  if (nrow(exon_found) == 0) return(list(coding_pos = NA, exon = NA))
  coding_pos <- exon_found$cumulative_length + exon_found$pos_within_exon
  return(list(coding_pos = coding_pos, exon = exon_found$exon))
}

# Plot function for variants, exons, and pathogenicity
variant_plot <- function(data_all, prefix = "variant_plot") {
  if (!"start.x" %in% colnames(data_all)) stop("🚨 Missing 'start.x' column.")
  if (!"end.x" %in% colnames(data_all)) stop("🚨 Missing 'end.x' column.")
  if (!"PHRED" %in% colnames(data_all)) stop("🚨 Missing 'PHRED' column.")
  if (!"pathogenic" %in% colnames(data_all)) stop("🚨 Missing 'pathogenic' column.")  
  
  # Compute c. positions and map pathogenic values
  data_all <- data_all %>%
    mutate(start = as.numeric(start.x),
           end = as.numeric(end.x),
           pathogenic = factor(case_when(
             pathogenic == -1 ~ "Benign",
             pathogenic == 0 ~ "Uncertain",
             pathogenic == 1 ~ "Pathogenic",
             TRUE ~ "Unknown"  # 🛡️ Safety fallback
           ), levels = c("Benign", "Uncertain", "Pathogenic", "Unknown"))) %>%
    rowwise() %>%
    mutate(cdot_info = list(find_cdot(start, end))) %>%
    ungroup() %>%  
    mutate(c_position = map_dbl(cdot_info, ~ .x$coding_pos),
           exon = factor(map_chr(cdot_info, ~ as.character(.x$exon)))) %>%
    dplyr::select(-cdot_info)  # ✅ Now works without error
  
  # Plot with exon boundaries
  p <- ggplot(data_all, aes(x = c_position, y = PHRED, color = pathogenic)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = c("Benign" = "green", "Uncertain" = "blue", "Pathogenic" = "red", "Unknown" = "gray")) +
    geom_vline(xintercept = gene_info$cumulative_length + 1, color = "black", linetype = "dotted") +
    theme_minimal(base_size = 14) +
    labs(title = "Variant Pathogenicity Across Coding Sequence",
         x = "c. Position (Coding DNA Position)",
         y = "PHRED Score",
         color = "Pathogenicity") +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"))
  
  # Save plot
  ggsave(paste0(prefix, "_pathogenicity_plot.pdf"), p, width = 10, height = 6)
  cat("✅ Plot saved as:", paste0(prefix, "_pathogenicity_plot.pdf"), "\n")
}

# Usage Example
variant_plot(d_all, prefix = "final_variant_plot")
