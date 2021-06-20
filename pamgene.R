#!/usr/local/bin/RScript

suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(broom)))
suppressWarnings(suppressMessages(library(MESS))) # for auc function
suppressWarnings(suppressMessages(library(sva))) # for combat function
suppressWarnings(suppressMessages(library(matrixStats)))
suppressWarnings(suppressMessages(library(modelr)))
suppressWarnings(suppressMessages(library(gtools)))
suppressWarnings(suppressMessages(library(parallel)))
suppressWarnings(suppressMessages(library(tidyverse)))

# functions -------------
# maybe better scaling function
#redist.fun <- function(x){(x-min(x))/diff(range(x))}

# maybe better scaling function
redist.fun <- function(x, no_perm){(x-0)/log10(no_perm)}

sum_peptide_score = function(x) {
  as.data.frame(x) %>% 
    dplyr::group_by(id) %>% 
    dplyr::summarise(pn = sum(score))
}

calculate_peptide_scores = function(pep_identifiers, value_matrix, info_df, ctrl_cols, exp_cols, center = TRUE, scale = TRUE) {
  df_scaled = scale(value_matrix, center = center, scale = scale)
  df_means = data.frame(id = pep_identifiers, 
                        A1 = rowMeans(df_scaled[,ctrl_cols]), 
                        A2 = rowMeans(df_scaled[,exp_cols]), 
                        A1s = rowSds(df_scaled[,ctrl_cols], na.rm = T), 
                        A2s = rowSds(df_scaled[,exp_cols], na.rm = T))
  df_means = df_means %>% 
    dplyr::mutate(pep_score = (A2 - A1)/sqrt(A1s^2 + A2s^2))
  phos_pred = merge(info_df, df_means, by = "id")
  # optional weighing by log prediction score
  # phos_pred = phos_pred %>% dplyr::mutate(pep_kinase_score = pep_score * log_pred_score)
}

estimate_top_kinases = function(phos_pred_df, kinexus_score, min_kinases = 5, max_kinases = 50) {
  diff_score = lapply(min_kinases:max_kinases, function(x) {
    phos_pred_filtered = phos_pred_df %>% 
      dplyr::filter(kinexus_pred_score_v2 >= kinexus_score) %>% 
      dplyr::filter(kinase_rank <= x)  %>% 
      dplyr::group_by(kinases_id) %>% 
      dplyr::summarise(diff_score = sum(pep_score)/n())
    phos_pred_filtered = colMeans(phos_pred_filtered[,2]) %>% 
      unlist() %>% 
      unname()
  })
  top_kinases = (min_kinases:max_kinases)[which.max(diff_score)]
}

predict_kinases = function(phos_pred_df, top_kinases) {
  kinase_list = phos_pred_df %>% 
    dplyr::filter(kinase_rank <= top_kinases) %>%
    dplyr::group_by(kinases_id) %>% 
    dplyr::mutate(n = n()) %>%
    dplyr::mutate(pep_norm = pep_score/n)
}

calculate_permutation_score = function(long_perm_df, kinase_summary, no_perm) {
  merge(long_perm_df, kinase_summary, by = "kinases_id") %>% 
    dplyr::mutate(exceeds = ifelse((abs(sum_score.x) - abs(sum_score.y)) > 0, 1, 0)) %>% 
    dplyr::group_by(kinases_id) %>% 
    dplyr::summarise(m_stat = sum(exceeds)) %>%
    rowwise %>% 
    dplyr::mutate(q = -log10(max(m_stat/no_perm, 1/no_perm)))
}

# parse commandline arguments -------------
option_list = list(
  make_option(c("-i", "--file_in"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--file_out"), type="character", default="kinase_scores.txt", 
              help="output file name [default=%default]", metavar="character"),
  make_option(c("-c", "--ctrl"), type="character", default="1,2,3,4", 
              help="comma delimited numbers of control columns in input file [default=%default]", metavar="character"),
  make_option(c("-e", "--exp"), type="character", default="5,6,7,8", 
              help="comma delimited numbers of experimental columns in input file [default=%default]", metavar="character"),
  make_option(c("--top_kinases"), type="character", default=NULL, 
              help="top n kinases identified for peptide", metavar="character"),
  make_option(c("--kinexus_score"), type="numeric", default=300, 
              help="min kinexus score for kinase-peptide combination [default=%default]"),
  make_option(c("--no_perm"), type="numeric", default=1000, 
              help="number of permutations for sensitivity and specificity calculations [default=%default]"),
  make_option(c("-b", "--batch_correction"), action = "store_true", default  = FALSE,
              help="apply chip batch correction by sva::ComBat [default=%default]")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# analysis options  -------------
center.all = T
scale.all = T
batch_corr = opt$batch_correction
no_perm = opt$no_perm
top_kinases = opt$top_kinases
kinexus_score = opt$kinexus_score

# load db -------------
setwd("~/OneDrive/programming/pamgene/")

# how to generate db rds files from csvs
# peptide_db = read.csv("dbs/peptide_db.csv",
#                       sep = ";",
#                       stringsAsFactors=FALSE) %>%
#   dplyr::select(uniprot_id, uniprot_name, identity, id, sequence, phosphosite) %>%
#   distinct()
# 
# phosphonet_db = read_delim("dbs/phosphonet_db.csv",
#                             delim = ";",
#                             escape_double = FALSE,
#                             trim_ws = TRUE)
# kinase_name_id = unique(phosphonet_db[,c(7,6)]) %>%
#   setNames(c("id", "name"))
# 
# peptide_phosphonet = merge(phosphonet_db, peptide_db, by.x = c("substrate", "position"), by.y = c("uniprot_id", "phosphosite"))
# saveRDS(kinase_name_id, "dbs/kinase_name_id.rds")
# saveRDS(peptide_phosphonet, "dbs/peptide_phosphonet.rds")

peptide_phosphonet = readRDS("dbs/peptide_phosphonet.rds") %>%  
  filter(kinexus_pred_score_v2 >= kinexus_score) %>%
  group_by(kinases_id, id) %>% 
  arrange(kinases_id, kinase_rank) %>%
  distinct(kinases_id, id, .keep_all = T) %>% 
  group_by(kinases_id) %>%
  dplyr::select(c(kinase_rank, kinase_name, position, kinexus_pred_score_v2, kinases_id, id))

kinase_name_id = readRDS("dbs/kinase_name_id.rds")

# experimental design  -------------
# raw_data_file = "raw_data_stk.csv"
# output_file = "kinase_scores.txt"
# ctrl_cols = c(1,2,3,4)
# exp_cols = c(5,6,7,8)
raw_data_file = opt$file_in
output_file = opt$file_out
ctrl_cols = str_split(opt$ctrl, ",")[[1]] %>% as.numeric()
exp_cols = str_split(opt$exp, ",")[[1]] %>% as.numeric()
no_samples = length(c(ctrl_cols, exp_cols))

# load data -------------
raw_data = read_delim(raw_data_file, ";", escape_double = FALSE, trim_ws = TRUE, col_types = cols())
raw_long = melt(raw_data) %>%
  separate(col = variable, into = c("condition", "chip", "time"), sep = "_") %>%
  dplyr::select(-UniprotAccession) %>%
  setNames(c("id", "condition", "chip", "time", "measurement")) %>%
  dplyr::mutate(time = as.numeric(time))

# calculate AUC -------------
# replaced in favour of linear model
# kinase_auc = raw_long %>% 
#   dplyr::group_by(id, condition, chip) %>% 
#   do(auc = auc(.$time, .$measurement)) %>%
#   tidyr::unnest(auc) %>%
#   dcast(id ~ condition + chip, value.var = "auc") # replace with pivot wider

kinase_auc = raw_long %>% 
  dplyr::filter(measurement < max(measurement)) %>% 
  group_by(id, condition, chip) %>% 
  do(lm = lm(.$measurement ~ .$time)) %>% 
  tidy(lm) %>%
  dplyr::filter(term == ".$time") %>%
  dplyr::select(id, condition, chip, estimate) %>%
  dcast(id ~ condition + chip, value.var = "estimate") %>%
  drop_na()

# ComBat batch effect normalization -----------
if(batch_corr == TRUE) {
  cat("Correcting for chip batch effects...\n")
  cb_metadata = raw_long %>% dplyr::select(condition, chip) %>% 
    distinct() %>%
    dplyr::mutate(sample = paste0(condition, "_", chip)) %>% 
    dplyr::select(-condition) %>%
    dplyr::select(sample, chip)
  
  cb_data = as.matrix(kinase_auc[,-1])
  rownames(cb_data) = kinase_auc[,1]
  
  cb_corr_model = model.matrix(~1, data = cb_metadata)
  cb_corr_counts = ComBat(dat=cb_data,
                          batch=cb_metadata$chip,
                          mod=cb_corr_model,
                          par.prior=TRUE,
                          prior.plot=FALSE)
  kinase_auc = cb_corr_counts %>% as_tibble(rownames = "id")
}

# TODO calculate peptide statistics -------------
tt = apply(kinase_auc[,2:(no_samples+1)],
            MARGIN = 1,
            function (x) {t.test(x[ctrl_cols], x[exp_cols], paired = F)})
pvals = unname(sapply(tt, function(x) {unlist(x[3])}))

# predict kinases -------------
peptide_scores = calculate_peptide_scores(kinase_auc[,1],
                                          kinase_auc[,2:(no_samples+1)],
                                          peptide_phosphonet, 
                                          ctrl_cols,
                                          exp_cols,
                                          center = center.all, 
                                          scale = scale.all)

if(is.null(top_kinases)) {
  cat("Determining optimal number of kinases to include...\n")
  top_kinases = estimate_top_kinases(peptide_scores,
                                     kinexus_score,
                                     min_kinases = 5, 
                                     max_kinases = 15)
  cat(paste0(top_kinases, " top kinases selected.\n"))
}

kinase_pred_summary = predict_kinases(peptide_scores, top_kinases) %>% 
  group_by(kinases_id) %>% 
  summarise(sum_score = sum(pep_norm),
            mean_score = mean(pep_norm),
            median_score = median(pep_norm),
            sd_score = sd(pep_norm), 
            no_peptides = mean(n))

cat(paste0("Number of permutations - ", no_perm, "\n"))
# specificity score - peptide permutation -------------
cat("Calculating specificity score - permuting peptides...\n")
peptide_perm_model = suppressWarnings(modelr::permute(kinase_auc, no_perm, columns = id))
peptide_perm = map_dfr(peptide_perm_model$perm, function(x) {
  y = as.data.frame(x)
  w = calculate_peptide_scores(y[,1], 
                               y[,2:(no_samples+1)], 
                               peptide_phosphonet, 
                               ctrl_cols, 
                               exp_cols, 
                               center = center.all, 
                               scale = scale.all)
  z = predict_kinases(w, top_kinases) %>% 
    group_by(kinases_id) %>% 
    summarise(sum_score = sum(pep_norm),
              mean_score = mean(pep_norm),
              sd_score = sd(pep_norm))
})
specificity_df = calculate_permutation_score(peptide_perm, kinase_pred_summary, no_perm) 

# selectivity score - sample_permutation -------------
cat("Calculating selectivity score - permuting samples...\n")
sample_perm = mclapply(1:no_perm, function(x) {
  w = calculate_peptide_scores(kinase_auc[,1],
                               kinase_auc[,2:(no_samples+1)][,gtools::permute(1:no_samples)],
                               peptide_phosphonet, 
                               ctrl_cols,
                               exp_cols,
                               center = center.all,
                               scale = scale.all)
  
    z = predict_kinases(w, top_kinases) %>% 
      group_by(kinases_id) %>% 
      summarise(sum_score = sum(pep_norm),
                mean_score = mean(pep_norm),
                sd_score = sd(pep_norm))
}) %>% bind_rows()

selectivity_df = calculate_permutation_score(sample_perm, kinase_pred_summary, no_perm)

# final kinase scores -------------
kinase_scores = cbind.data.frame(kinase_pred_summary,
                                 specificity = redist.fun(specificity_df$q, no_perm), 
                                 selectivity = redist.fun(selectivity_df$q, no_perm)) %>%
  dplyr::mutate(total = specificity + selectivity) %>%
  dplyr::left_join(kinase_name_id, by = c("kinases_id" = "id")) %>%
  dplyr::select(kinases_id, name, specificity, selectivity, total, everything()) %>%
  dplyr::arrange(-total, -abs(sum_score)) %>%
  drop_na

# export data -------------
write_delim(kinase_scores, 
            path = output_file, 
            delim = "\t")

cat("Done.\n")
