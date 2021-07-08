#library(pheatmap)
#library(tidyverse)

# basic visualisation -------------

# heatmap - compare overall peptide signals
pheatmap(scale(kinase_auc[,2:9], center = T, scale = F), 
         cluster_cols = T, 
         cluster_rows = T, 
         scale = "row")

# boxplot - compare overall peptide signals
kinase_auc %>% 
  ggplot(aes(x = condition, y = estimate)) + 
  geom_boxplot() + 
  facet_wrap(vars(id), ncol = 10, scales = "free_y")

# graph data from pathway analysis
# filtered StringDB pathway data is accepted - GO, Reactome, KEGG, Uniprot keywords
pathway_analysis_file = "~/OneDrive/_Ruaslab/01_projects/13_lmcd1_duarte/pamgene/pamgene_8_5/pamgene_summary_8_5_filtered.xlsx"
sheet_no = gdata::sheetCount(pathway_analysis_file)

pathway_analysis = lapply(1:sheet_no, function(x) {
  df = read.xlsx(pathway_analysis_file, sheetIndex = x)
})

# kinases graph
upstream_analysis[[1]] %>% arrange(t) %>% mutate(name = factor(name, name)) %>% 
  ggplot(aes(x = name, y = t, fill = specificity)) + geom_bar(stat="identity") + coord_flip() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + theme_bw()
# TODO add ggsave

# kinases graph - walk through sheets that contain pathway analysis
walk(3:5, function(x) {
  upstream_analysis[[x]] %>% arrange(log10pval) %>% mutate(pathway = factor(pathway, pathway)) %>% 
    ggplot(aes(x = pathway, y = log10pval)) + geom_bar(stat="identity") + coord_flip() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + theme_void()
  # TODO add ggsave
})
