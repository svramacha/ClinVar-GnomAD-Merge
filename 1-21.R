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
  #aucs <- numeric(n_imputations) # Store AUC imputation

  
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
              #aucs =AUCs, 
              #roc_data = roc_data, 
              plot_data = plot_data))
  }











### STOPPED HERE









#auc_imputation <- function(d_obs, d_unc, annot, p_train=.8, n_imp=200){
  
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
    
    return(cutoff - odds_bound / odds_path)
  }
  
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
                              "benign", cutlabel[j], result$root, rng))
    }       
  }
  
  rslts <- as.data.frame(rslts, stringsAsFactors=F)
  names(rslts) <- c("annot", "int", "logor", "se_int", "se_logor", "pval", "var_type", 
                    "acmg_cat", "annot_value", "min_obs", "max_obs")
  
  tbl <- rslts %>% tidyr::pivot_wider(names_from=c("var_type","acmg_cat"),
                                      values_from="annot_value",
                                      names_glue = "{var_type}_{acmg_cat}") %>%
    dplyr::select(annot, int, logor, se_int, se_logor, pval, 
                  starts_with("pathogenic"), starts_with("benign"),
                  min_obs, max_obs)
  
  
  write.table(file = file, tbl, quote=F, row.names=F, col.names=T, sep="\t")
  cat("Wrote threshold for different ACMG support categories to", file, "\n")
  
  return(tbl)
}


    
    
    
