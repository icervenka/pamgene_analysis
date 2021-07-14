#!/usr/local/bin/RScript

# imports -------------
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(broom)))
suppressWarnings(suppressMessages(library(MESS))) # for auc function
suppressWarnings(suppressMessages(library(matrixStats)))
suppressWarnings(suppressMessages(library(modelr)))
suppressWarnings(suppressMessages(library(gtools)))
suppressWarnings(suppressMessages(library(purrr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(tibble)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(dplyr)))

# functions -------------

# TODO maybe better scaling function
#' Scale kinase selectivity and specificity measures to 0-1 range
#' based on number of permutations
#'
#' @param x measure to scale
#' @param no_perm number of permutations used to calculate x
#'
#' @return x scaled do (0,1) range
#'
#' @examples
redist.fun = function(x, no_perm){(x-0)/log10(no_perm)}

#' Calculates peptide scores based on mean and sd of control and experimental groups
#'
#' @param data_df data.frame of peptide kinetics measurements where first column
#' is the peptide id. 
#' @param metadata named list created by 'parse_sample_metadata' function, has
#' to contain $ctrl and $exp names containing column positions of control and
#' exprimental samples
#' @param info_df data.frame summarising the pamgene peptide chip data and 
#' phosphonet kinase data
#' @param ctrl_cols integer vector of control column position in peptide
#' kinetics matrix 
#' @param exp_cols integer vector of experimental column position in peptide
#' kinetics matrix 
#' @param center logical value - center the matrix of peptide values
#' @param scale logical value - scale the matrix of peptide values
#'
#' @return data.frame with calculated score for each peptide on pamgene chip
#' based on following formula:
#' (mean(exp_group) - mean(ctrl_group))/sqrt(sd(ctrl_group)^2 + sd(exp_group)^2)
#'
#' @examples
#' TODO move merging into another function
#' TODO provide data as single variable and create matrix with rownames instead
calculate_peptide_scores = function(data_df, metadata, info_df, 
                                    center = TRUE, scale = TRUE) {
  # convert data.frame to matrix with rownames for scaling
  data_matrix = as.matrix(data_df[-1])

    rownames(data_matrix) = data_df %>% dplyr::pull(id)
  df_scaled = scale(data_matrix, center = center, scale = scale)
  
  # peptide score calculation of scaled and centered data
  df_means = data.frame(id = rownames(data_matrix),
                        A1 = rowMeans(df_scaled[,metadata$ctrl]),
                        A2 = rowMeans(df_scaled[,metadata$exp]),
                        A1s = matrixStats::rowSds(df_scaled[,metadata$ctrl], na.rm = T),
                        A2s = matrixStats::rowSds(df_scaled[,metadata$exp], na.rm = T))
  df_means = df_means %>%
    dplyr::mutate(pep_score = (A2 - A1)/sqrt(A1s^2 + A2s^2))
  phos_pred = merge(info_df, df_means, by = "id")
  # optional weighing by log prediction score
  # phos_pred = phos_pred %>% dplyr::mutate(pep_kinase_score = pep_score * log_pred_score)
}

#' Parses group information for samples from commandline arguments. Group info
#' can be specified either as comma-separated indices of alphabetically 
#' arranged <group_name>_<chip_id> denominators or a single <group_name> string.
#' For the calculations of column indices from <group_name> first column containing
#' peptide ids is ommited.
#'
#' @param opt commandline arguments object
#' @param peptide_reaction_estimate data.frame containing peptide reaction 
#' velocity estimates with first column being id and subsequent columns belonging
#' each to one sample
#'
#' @return named list containing three items:
#' - .$ctrl column position of control samples
#' - .$exp column positions of experimental samples
#' - .$n total number of samples
#'
#' @examples
parse_sample_metadata = function(opt, peptide_reaction_estimate) {
  # process control and experimental groups commandline args
  ctrl_cols = stringr::str_split(opt$ctrl, ",")[[1]] 
  exp_cols = stringr::str_split(opt$exp, ",")[[1]]
  
  # transform to vector of column positions (excluding the id column)
  if(length(c(ctrl_cols, exp_cols)) == 2) {
    ctrl_cols = which(grepl(ctrl_cols, names(peptide_reaction_estimate[-1])))
    exp_cols = which(grepl(exp_cols, names(peptide_reaction_estimate[-1])))
  } else {
    ctrl_cols = ctrl_cols %>% as.numeric()
    exp_cols = exp_cols %>% as.numeric()
  }
  
  return(list("ctrl" = ctrl_cols, 
              "exp" = exp_cols, 
              "n" = length(c(ctrl_cols, exp_cols))))
}

#' Estimates the number of top kinases to use in kinase prediction. Phosphonet 
#' website  reports 50 kinases as being able to phosphorylate specific peptide, 
#' using only subset of top N kinases can give better results. This function is 
#' invoked if the number of top kinases is not specified on commandline. 
#' It works by iteratively increasing the number of included kinases, 
#' calculating the mean of average peptide scores for each set of kinases, 
#' and choosing the group with maximum value
#'
#' @param phos_pred_df data.frame of peptide scores and kinases reported for 
#' each pamgene peptide on the chip. 
#' @param kinexus_score integer - threshold kinexus v2 score the kinase has to 
#' pass to be considered. Low score kinases are filtered out before selecting 
#' top N kinases to include
#' @param min_kinases integer - minumum number of kinases to consider. 
#' The default value is 5
#' @param max_kinases integer - minumum number of kinases to consider. 
#' The default value is 50
#'
#' @return integer - number of top N kinases to include in prediction
#'
#' @examples
estimate_top_kinases = function(phos_pred_df, kinexus_score, 
                                min_kinases = 5, max_kinases = 50) {
  diff_score = lapply(min_kinases:max_kinases, function(x) {
    phos_pred_filtered = phos_pred_df %>%
      dplyr::filter(kinexus_score_v2 >= kinexus_score) %>%
      dplyr::filter(kinase_rank <= x)  %>%
      dplyr::group_by(kinase_id) %>%
      dplyr::summarise(diff_score = sum(pep_score)/n())
    
    # calculate average of diff_score
    phos_pred_filtered = colMeans(phos_pred_filtered[,2]) %>%
      unlist() %>%
      unname()
  })
  top_kinases = (min_kinases:max_kinases)[which.max(diff_score)]
}

#' Predicts kinases that phosphorylate specific pamgene peptide by dividing
#' the peptide score for each peptide by factor corresponding to the total number
#' of times the kinase appears among top N kinases.
#'
#' @param phos_pred_df peptide scores for pamgene peptides 
#' @param top_kinases top N kinases for each peptide to include in the 
#' calculations. Can be estimated from data by 'estimate_top_kinses' function or 
#' specified directly.
#'
#' @return
#'
#' @examples
predict_kinases = function(phos_pred_df, top_kinases) {
  kinase_list = phos_pred_df %>%
    dplyr::filter(kinase_rank <= top_kinases) %>%
    dplyr::group_by(kinase_id) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::mutate(pep_norm = pep_score/n) %>%
    dplyr::group_by(kinase_id) %>%
    dplyr::summarise(sum_score = sum(pep_norm),
                     mean_score = mean(pep_norm),
                     sd_score = sd(pep_norm))
}

#' Calculate final metrics aggregating calculation on permutated data (row or column).
#' First the number of times the peptide scores from all permutations exceeds the 
#' peptide score calculated on the original data set is reported as m_stat. Next,
#' q statistics expressed as -log10(max(m_stat/no_perm, 1/no_perm) is reported
#'
#' @param long_perm_df data.frame summarising the results from permutation calculations
#' @param kinase_summary data.frame of predicted kinases for phosphosites
#' @param no_perm integer, number of permutations used to calculate the statistics
#'
#' @return data.frame containing following columns:
#' - kinase_id
#' - m_stat
#' - q
#'
#' @examples
calculate_permutation_score = function(long_perm_df, kinase_summary, no_perm) {
  merge(long_perm_df, kinase_summary, by = "kinase_id") %>%
    dplyr::mutate(exceeds = ifelse((abs(sum_score.x) - abs(sum_score.y)) > 0, 1, 0)) %>%
    dplyr::group_by(kinase_id) %>%
    dplyr::summarise(m_stat = sum(exceeds), .groups = "drop") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(q = -log10(max(m_stat/no_perm, 1/no_perm)))
}

# parse commandline arguments -------------
option_list = list(
  make_option(c("-i", "--file_in"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--dir_out"), type="character", default=".",
              help="output directory [default=current dir]", metavar="character"),
  make_option(c("--pamgene"), type="character", default="pamgene.csv",
              help="list of peptides included in pamgene chip, see Readme.md 
              for formatting details [default=%default]", metavar="character"),
  make_option(c("--phosphonet"), type="character", default="phosphonet.csv",
              help="kinase prediction file for peptides included in pamgene chip.
              Can be generated using phosphonet scraper. See Readme.md for more 
              details. [default=%default]", metavar="character"),
  make_option(c("-c", "--ctrl"), type="character", default="1,2,3,4",
              help="comma delimited numbers of control columns in input file or a 
              single string indicating group name [default=%default]", 
              metavar="character"),
  make_option(c("-e", "--exp"), type="character", default="5,6,7,8",
              help="comma delimited numbers of experimental columns in input file 
              or asingle string indicating group name [default=%default]", 
              metavar="character"),
  make_option(c("--top_kinases"), type="character", default=NULL,
              help="top n kinases identified for peptide", metavar="character"),
  make_option(c("--kinexus_score"), type="numeric", default=300,
              help="min kinexus score for kinase-peptide combination 
              [default=%default]"),
  make_option(c("--no_perm"), type="numeric", default=1000,
              help="number of permutations for sensitivity and specificity 
              calculations [default=%default]"),
  make_option(c("-b", "--batch_correction"), action = "store_true", default  = FALSE,
              help="apply chip batch correction by sva::ComBat [default=%default]")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# analysis options and arguments  -------------
# center and scaling of the peptide values
# can be modified if centering and scaling of measurements is not desireable
# TODO optionally move to commandline parameters
center.all = T
scale.all = T

# assign command line args to variables
batch_corr = opt$batch_correction
no_perm = opt$no_perm
top_kinases = opt$top_kinases
kinexus_score = opt$kinexus_score
raw_data_file = opt$file_in
output_dir = opt$dir_out

# generate combined peptide info -------------
# read info on pamgene peptides
tryCatch({
  pamgene = read.csv(opt$pamgene,
                     stringsAsFactors = F) %>%
    dplyr::select(id,
                  uniprot_id, 
                  phosphosite) %>%
    dplyr::distinct()
}, error = function(e) {
  stop("Pamgene peptide file has incorrect format. See Readme.md for more details.")
})

# read info on predicted kinases from phosphonet.ca
tryCatch({
  phosphonet = read.csv(opt$phosphonet,
                        stringsAsFactors = F) %>%
    dplyr::select(substrate,
                  site,
                  kinase_rank,
                  kinase_name,
                  kinase_id,
                  kinexus_score_v2)
}, error = function(e) {
  stop("Phosphonet peptide file has incorrect format. See Readme.md for more details.")
})

# extract data.frame mapping of kinase ids to kinase names
kinase_name_id = phosphonet %>%
  dplyr::select(kinase_id, kinase_name) %>%
  unique()

# merge pamgene and phosphonet data to single 'metadata' data.frame
peptide_phosphonet = merge(phosphonet, 
                           pamgene, 
                           by.x = c("substrate", "site"), 
                           by.y = c("uniprot_id", "phosphosite")) %>%
  dplyr::filter(kinexus_score_v2 >= kinexus_score) %>%
  dplyr::group_by(kinase_id, id) %>%
  dplyr::arrange(kinase_id, kinase_rank) %>%
  dplyr::distinct(kinase_id, id, .keep_all = T) %>%
  dplyr::group_by(kinase_id) %>%
  dplyr::select(c(kinase_rank, 
                  kinase_name, 
                  site, 
                  kinexus_score_v2, 
                  kinase_id, 
                  id))

# load data -------------
# raw_data = readr::read_delim(raw_data_file, 
#                              ";", 
#                              escape_double = FALSE, 
#                              trim_ws = TRUE, 
#                              col_types = readr::cols())

raw_data = read.csv(raw_data_file,
                    stringsAsFactors = F)

# transform if not in long format already
if(dim(raw_data)[2] != 6) {
  raw_long = raw_data %>%
    tidyr::pivot_longer(-c(ID, UniprotAccession), 
                        names_to = "variable", 
                        values_to = "value") %>%
    tidyr::separate(col = variable, 
                    into = c("group", "chip_id", "time"), 
                    sep = "_") %>%
    dplyr::mutate(time = as.numeric(time))
} else if(names(raw_data) == c("ID", "UniprotAccession",  "group", "chip_id", "time", "value")) {
  raw_long = raw_data
} else {
  stop("Format of input data not recognized. Please reformat your data and try again.")
}

# using linear model fitting for the lack of better measurement
# such as reaction velocity
# removes saturated data points before fitting the model
peptide_reaction_estimate = raw_long %>%
  dplyr::filter(value < max(value)) %>%
  dplyr::group_by(ID, group, chip_id) %>%
  do(broom::tidy(lm(.$value ~ .$time))) %>%
  dplyr::filter(term == ".$time") %>%
  dplyr::select(id = ID, group, chip_id, estimate) %>%
  tidyr::pivot_wider(id_cols = "id", 
                     names_from = c("group", "chip_id"), 
                     values_from = "estimate") %>%
  tidyr::drop_na()

# ComBat batch effect normalization -----------
if(batch_corr == TRUE) {
  cat("Correcting for chip batch effects...\n")
  suppressWarnings(suppressMessages(library(sva))) # for combat function
  cb_metadata = raw_long %>% 
    dplyr::select(condition, chip) %>%
    dplyr::distinct() %>%
    dplyr::mutate(sample = paste0(condition, "_", chip)) %>%
    dplyr::select(-condition) %>%
    dplyr::select(sample, chip)

  cb_data = as.matrix(peptide_reaction_estimate[-1])
  rownames(cb_data) = peptide_reaction_estimate[1]

  cb_corr_model = model.matrix(~1, data = cb_metadata)
  cb_corr_counts = ComBat(dat=cb_data,
                          batch=cb_metadata$chip,
                          mod=cb_corr_model,
                          par.prior=TRUE,
                          prior.plot=FALSE)
  peptide_reaction_estimate = cb_corr_counts %>% as_tibble(rownames = "id")
}

# parse sample metadata from commandline -------------
# currently needs peptide_reaction_estimate to map column indices to groups
# this info can be inferred from raw_data, I can have a look at it later
meta = parse_sample_metadata(opt, peptide_reaction_estimate)

# predict kinases -------------
# calculate peptide scores
peptide_scores = calculate_peptide_scores(peptide_reaction_estimate,
                                          meta,
                                          peptide_phosphonet,
                                          center = center.all,
                                          scale = scale.all)

# if top_kinases are not specified on commandline, predict from provided data
# pamgene recommends not using more than 15, but the value can be overridden
if(is.null(top_kinases)) {
  cat("Determining optimal number of kinases to include...\n")
  top_kinases = estimate_top_kinases(peptide_scores,
                                     kinexus_score,
                                     min_kinases = 5,
                                     max_kinases = 15)
  cat(paste0(top_kinases, " top kinases selected.\n"))
}

# use peptide scores to predict most likely reponsible kinases
kinase_pred_summary = predict_kinases(peptide_scores, top_kinases)

# specificity score - peptide permutation -------------
cat("Calculating specificity score - permuting peptides...\n")
# progress bar
pbp = txtProgressBar(min = 0,
                     max = no_perm,
                     style = 1,
                     width = 70,
                     char = "=")

# create data.frame of row permutations to feed into the calculations
peptide_perm_model = suppressWarnings(modelr::permute(peptide_reaction_estimate, 
                                                      no_perm, 
                                                      columns = id))

# iterate through permutations, calculating peptide scores together
# together with kinase prediction for each 
# return the combined data.frame of kinase prediction for all permutations
peptide_perm = purrr::map_dfr(1:no_perm, function(x) {
  setTxtProgressBar(pbp, x)
  y = as.data.frame(peptide_perm_model[[1]][[x]])
  w = calculate_peptide_scores(y,
                               meta,
                               peptide_phosphonet,
                               center = center.all,
                               scale = scale.all)
  z = predict_kinases(w, top_kinases)
})
close(pbp)

# merge together prediction of kinases by counting how many times
# the score from original data is higher then the one from permutations
specificity_df = calculate_permutation_score(peptide_perm, 
                                             kinase_pred_summary, 
                                             no_perm) %>%
  dplyr::ungroup()

# selectivity score - sample_permutation -------------
cat("Calculating selectivity score - permuting samples...\n")
# progress bar
pbs = txtProgressBar(min = 0,
                     max = no_perm,
                     style = 1,
                     width = 70,
                     char = "=")
# Here the data.frame of permutations is not calculated beforehand, the column
# positions for permutations are passed directly into calculate_peptide_scores 
# function return the combined data.frame of kinase prediction for all permutations
sample_perm = purrr::map_dfr(1:no_perm, function(x) {
  setTxtProgressBar(pbs, x)
  # permutation of column numeric indices with id column added back 
  column_perm = data.frame(peptide_reaction_estimate[1], 
                           peptide_reaction_estimate[-1][,gtools::permute(1:meta$n)])
  w = calculate_peptide_scores(column_perm,
                               meta,
                               peptide_phosphonet,
                               center = center.all,
                               scale = scale.all)

  z = predict_kinases(w, top_kinases)
})
close(pbs)

# merge together prediction of kinases by counting how many times
# the score from original data is higher then the one from permutations
selectivity_df = calculate_permutation_score(sample_perm, 
                                             kinase_pred_summary, no_perm) %>%
  dplyr::ungroup()

# final kinase scores -------------
# combine selectivity and specificity data.frames created by row and column 
# permutations rearrange the data.frame columns and order rows for export
kinase_scores = cbind.data.frame(kinase_pred_summary,
                                 # these measures are always in range (0, log10(no_perm)) 
                                 # and need to be scaled to (0,1)
                                 specificity = redist.fun(specificity_df$q, no_perm),
                                 selectivity = redist.fun(selectivity_df$q, no_perm)) %>%
  dplyr::mutate(total = specificity + selectivity) %>%
  dplyr::left_join(kinase_name_id, by = "kinase_id") %>%
  dplyr::select(kinase_id, 
                kinase_name, 
                specificity, 
                selectivity, 
                total, 
                everything()) %>%
  dplyr::arrange(-total, -abs(sum_score)) %>%
  tidyr::drop_na()

# calculate peptide statistics -------------
# t.test between control and experimental groups, use broom::tidy to
# extract htest objects into the dataframe, merge with individual 
# sample-peptide measurements
# TODO add adjusted pvalue
peptide_stats = peptide_reaction_estimate %>%
  tibble::column_to_rownames(var = "id") %>%
  dplyr::group_by(id = rownames(.)) %>%
  do(broom::tidy(t.test(.[meta$ctrl], .[meta$exp]))) %>%
  dplyr::ungroup() %>%
  dplyr::select(id,
                mean_ctrl = estimate1,
                mean_exp = estimate2,
                statistic = statistic,
                pvalue = p.value) %>%
  dplyr::left_join(peptide_reaction_estimate, by = "id")

# export data -------------
cat("Exporting peptide scores.\n")
write.table(peptide_stats,
            file = paste0(output_dir, "/", "peptide_scores.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat("Exporting kinase scores.\n")
write.table(kinase_scores,
            file = paste0(output_dir, "/", "kinase_scores.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat("Done.\n")
