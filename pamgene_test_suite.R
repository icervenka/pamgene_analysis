peptide_num = sapply(kinase_name_id$name, function(x) {
  dim(peptide_phosphonet %>% filter(kinase_name == x))[1]
})
pep_nos2 = data.frame(name = names(pep_nos2), count = pep_nos2) %>% arrange(desc(pep_nos2))

peptide_num = peptide_phosphonet %>%
  dplyr::group_by(kinase_name) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::arrange(-count)



test_kinase = function(kinases_id, df, peptide_phosphonet, top_kinases, kinexus_score, cols, range_low, range_high, is_in_top_no) {
  test_stk = df
  #kinase_pep_ids = stk_pep_phosphonet %>% filter(`Kinases ID` == kinases_id, `Kinase Rank` < top_kinases, `Kinexus Predictor Version 2 Score` > kinexus_score) %>% dplyr::select(ID) %>% dplyr::distinct()
  kinase_pep_ids = peptide_phosphonet %>% filter(`Kinases ID` == kinases_id) %>% dplyr::select(ID) %>% dplyr::distinct()
  
  dim_all = dim(test_stk)
  dim_replace = dim(test_stk[test_stk$id %in% kinase_pep_ids$ID,cols])
  
  test_stk[,2:dim_all[2]] = sample(range_low, dim_all[1]*dim_all[2], replace = T)
  test_stk[test_stk$id %in% kinase_pep_ids$ID,cols] = sample(range_high, dim_replace[1]*dim_replace[2], replace = T)
  
  st_phos_pred = calculate_peptide_scores(test_stk[,1], test_stk[,2:(no_samples+1)],  peptide_phosphonet, ctrl_cols, exp_cols, center = center.all, scale = scale.all)
  st_kinase = predict_kinases(st_phos_pred, top_kinases)
  st_kinase_summary = st_kinase %>% group_by(id) %>% summarise(pn = sum(score))
  
  st_pep_perm = modelr::permute(test_stk, no_perm, columns = id)
  st_pep_perm = map(st_pep_perm$perm, function(x) {
    y = as.data.frame(x)
    w = calculate_peptide_scores(y[,1], y[,2:(no_samples+1)], peptide_phosphonet, ctrl_cols, exp_cols, center = center.all, scale = scale.all)
    z = predict_kinases(w, top_kinases) %>% group_by(id) %>% summarise(pn = sum(score))
  })
  st_pep_perm = bind_rows(st_pep_perm)
  
  specificity_df = calculate_permutation_score(st_pep_perm, st_kinase_summary, no_perm) %>% arrange(desc(q))
  
  ord = st_kinase_summary %>% arrange(desc(pn))
  pos_orig = match(kinases_id, ord$id)
  pos_orig_sum = pos_orig_sum + pos_orig
  print(ord[1,])
  print(peptide_phosphonet %>% filter(`Kinases ID` == ord[1,]$id))
  print(ord[pos_orig,])
  print(peptide_phosphonet %>% filter(`Kinases ID` == ord[pos_orig,]$id))
  test_top = specificity_df$id[1:is_in_top_no]
  pos = match(kinases_id, specificity_df$id)
  ret = ifelse(kinases_id %in% test_top, T, F)
  q = specificity_df[specificity_df$id == kinases_id, ]['q']
  print(paste(kinases_id, q, ret, pos_orig, pos, dim(peptide_phosphonet %>% filter(`Kinases ID` == kinases_id))[1]))
  print(specificity_df)
  print(pos_orig_sum)
  pos_orig_sum
}

test_top = sapply(kinase_name_id$id[1:100], function(x) {
  test_kinase(x, rls_cast_2, peptide_phosphonet, top_kinases, 300, c(6:9), c(1:5), c(20:40), 10)
})

rls_cast_2 = rls_cast
per = data.frame(permutations(152, 2)) %>% filter(X1 == 8 & X2 >= 5 & X2 <= 15)
maximize_pn_5 = lapply(1:dim(per)[1], function(x) {
  
  phosphonet_db = phosphonet_db %>% filter(`Kinase Rank` < per[x, 1], `Kinexus Predictor Version 2 Score` > kinexus_score)
  peptide_phosphonet = merge(phosphonet_db, stk_pep, by.x = "Substrate", by.y = "UniProt_ID")
  peptide_phosphonet = peptide_phosphonet %>% group_by(`Kinases ID`, ID) %>% arrange(`Kinases ID`, desc(`Kinexus Predictor Version 2 Score`)) %>%
    distinct(`Kinases ID`, ID, .keep_all = T) %>% group_by(`Kinases ID`) %>% top_n(per[x, 2], wt = ID) %>% dplyr::select(c(5,6,7,12))
  pos_orig = 0
  
  sum_positions = sapply(kinase_name_id$id, function(y) {
    kinase_pep_ids = peptide_phosphonet %>% filter(`Kinases ID` == y) %>% dplyr::select(`Kinases ID`, ID) %>% dplyr::distinct()
    
    test_stk = rls_cast_2
    dim_all = dim(test_stk)
    
    dim_replace = dim(test_stk[test_stk$id %in% kinase_pep_ids$ID,c(6:9)])
    if(dim_replace[1] > 0) {
      test_stk[,2:dim_all[2]] = sample(c(1:5), dim_all[1]*dim_all[2], replace = T)
      test_stk[test_stk$id %in% kinase_pep_ids$ID,c(6:9)] = sample(c(20:40), dim_replace[1]*dim_replace[2], replace = T)
      
      st_phos_pred = calculate_peptide_scores(test_stk[,1], test_stk[,2:(no_samples+1)],  peptide_phosphonet, ctrl_cols, exp_cols, center = center.all, scale = scale.all)
      st_kinase = predict_kinases(st_phos_pred, top_kinases)
      st_kinase_summary = st_kinase %>% group_by(id) %>% summarise(pn = sum(score))
      
      ord = st_kinase_summary %>% arrange(desc(pn))
      pos_orig = pos_orig + match(y, ord$id)
    }
    pos_orig
  })
  print(paste(per[x, 1], per[x, 2]))
  avg = data.frame(sun = sum_positions, kinases = dim(unique(phosphonet_db[,c(7,6)]))[1], X1 = per[x, 1], X2 = per[x, 2])
  avg
})

maximize_sum = bind_rows(maximize_pn_5)
maximize_sum = maximize_sum %>% group_by(X1, X2, kinases) %>% summarize(sum = sum(sun), kin = mean(kinases), zeros = length(which(sun == 0)), total_kin = n())
maximize_sum = maximize_sum %>% mutate(coverage = (total_kin-zeros)/total_kin, avg_pos = sum/kin) %>% arrange(avg_pos)


#Q13873 Q86Z02
kinase_pep_ids = peptide_phosphonet %>% filter(`Kinases ID` == "P78527", `Kinase Rank` <= 8, `Kinexus Predictor Version 2 Score` > 300) %>% dplyr::select(ID) %>% dplyr::distinct()
peptide_phosphonet %>% filter(`Kinases ID` == "P78527", `Kinase Rank` <= 8, `Kinexus Predictor Version 2 Score` < 300) %>% dplyr::select(ID) %>% dplyr::distinct()
peptide_phosphonet %>% filter(`Kinases ID` == "P45983", `Kinase Rank` <= 8, `Kinexus Predictor Version 2 Score` < 300) %>% dplyr::select(ID) %>% dplyr::distinct()

st_phos_pred %>% filter(ID == "FOXO3_25_37")
st_kinase %>% filter(id == "Q9UQ07")
st_kinase %>% filter(id == "P24941")

peptide_phosphonet %>% filter(`Kinases ID` == "Q9UQ07")
peptide_phosphonet %>% filter(`Kinases ID` == "P24941")

test_stk[,2:dim_data[2]] = sample(c(1:10), dim_data[1]*(dim_data[2]-1), replace = T)
test_stk[test_stk$id %in% kinase_pep_ids$ID,cols] = sample(c(20:50), length(cols)*dim_pep[1], replace = T)