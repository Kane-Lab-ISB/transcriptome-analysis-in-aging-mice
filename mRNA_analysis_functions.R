## Functions for subsequent analysis 

get_treat_contrast = function(tissue1, age1, sex1, tissue2, age2, sex2){
  contrast = paste0("group", tissue1, ".", age1, ".", sex1, " - ", "group", tissue2, ".", age2, ".", sex2)
  contrast.matrix = makeContrasts(contrasts = contrast, levels = mm)
  tmp = contrasts.fit(fit, contrast.matrix)
  tmp = treat(tmp, lfc = log2(1.5))
  tt = topTable(tmp, sort.by = "none", n = Inf, adjust.method = "BH")
  tt["gene"] = gene2de
  return(tt)
} 

get_volcano_plot = function(table){
  table_sig = table %>%
    filter(adj.P.Val < 0.05) %>%
    arrange(adj.P.Val)
  table_top10 = table_sig[1:10, ]
  table_sig_pos = table_sig[table_sig$logFC > 0, ]
  table_sig_neg = table_sig[table_sig$logFC < 0, ]
  plot = ggplot(data = table,
                aes(x = logFC,
                    y = -log10(P.Value))) +
    geom_point(size = 0.6, color = "grey95") +
    geom_point(data = table_sig_pos,
               aes(x = logFC,
                   y = -log10(P.Value)),
               color = "#39568CFF") +
    geom_point(data = table_sig_neg,
               aes(x = logFC,
                   y = -log10(P.Value)),
               color = "#DCE319FF") +
    geom_text(data = table_top10,
              aes(x = logFC,
                  y = -log10(P.Value),
                  label = table_top10$gene),
              vjust = -0.8,
              label.size = 2.5)
  return(plot)
}

heatmap_plt = function(tissue, table){
  de_table = table %>%
    filter(adj.P.Val < 0.05) 
  tissue_specific_samples = meta_dat_ids[meta_dat_ids$Tissue == tissue, ] %>%
    filter(Group == "old") 
  ## slice table and choose up to 25 genes by pvalue for heatmap
  overlap_genes_ms = norm_expr %>%
    as.data.frame() %>%
    dplyr::select(unlist(tissue_specific_samples$`Novogene ID`))
  overlap_genes_ms = overlap_genes_ms[rownames(overlap_genes_ms) %in% de_table$gene, ]
  overlap_genes_ms_scale = apply(overlap_genes_ms, 1, scale) %>%
    t()
  colnames(overlap_genes_ms_scale) = tissue_specific_samples$`Novogene ID`
  heatmap.2(overlap_genes_ms_scale, 
            trace = c("none"), 
            dendrogram = "column", 
            col = viridis(20), labRow = FALSE, 
            key.title = NA, key.xlab = NA, key.ylab = NA) 
}

enrich_analysis = function(gene_lst){
  dbs = c("KEGG_2019_Mouse") 
  enriched = enrichR::enrichr(gene_lst, dbs)
  enr_kegg = enriched$KEGG_2019_Mouse %>%
    visual_data_kegg() %>% 
    dplyr::select(Term, count, Adjusted.P.value, ratio)
  return(enr_kegg)
}

visual_data_kegg = function(table){
  table_out = table %>%
    arrange(Adjusted.P.value) %>%
    mutate(count = sapply(1:length(table$Overlap), 
                          function(x) {
                            as.numeric(strsplit(table$Overlap, "/")[[x]][1])
                          })) %>% 
    mutate(gene_all = sapply(1:length(table$Overlap), 
                             function(x) {
                               as.numeric(strsplit(table$Overlap, "/")[[x]][2])
                             })) %>%
    mutate(ratio = count/gene_all*100) %>%
    as.data.frame() 
  table_out = table_out[1:10, ] %>% 
    arrange(desc(Adjusted.P.value)) 
  return(table_out)
}

agecomp_plot = function(table1, table2){
  agecomp = merge(table1, table2, by = "gene")
  de_table1 = table1 %>%
    filter(adj.P.Val < 0.05)
  de_table2 = table2 %>%
    filter(adj.P.Val < 0.05)
  age_common = intersect(de_table1$gene, de_table2$gene)
  agecomp_sig_f = agecomp %>%
    filter(adj.P.Val.x < 0.05)
  agecomp_sig_m = agecomp %>%
    filter(adj.P.Val.y < 0.05)
  agecomp_com = agecomp[agecomp$gene %in% age_common, ]
  plot = ggplot(data = agecomp,
                aes(x = logFC.x,
                    y = logFC.y)) +
    geom_point(size = 0.6, color = "#E5E4E2") +
    geom_vline(xintercept = 0, lty = 3) +
    geom_hline(yintercept = 0, lty = 3) +
    labs(x = "log2 Fold Change - Females",
         y = "log2 Fold Change - Males") +
    geom_point(data = agecomp_sig_f, color = "#39568CFF", size = 3.5) +
    geom_point(data = agecomp_sig_m, color = "#DCE319FF", size = 3.5) +
    geom_point(data = agecomp_com, color = "#20A387FF", size = 5, shape = 17)
  return(plot)
}