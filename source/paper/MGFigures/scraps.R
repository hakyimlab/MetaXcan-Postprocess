build_data_by_pheno <- function() {
  d <- data.frame(pred_perf_r2=numeric(), 
                  zscore = numeric(),
                  pval = numeric(),
                  the_facet=character())
  for(phenoname in pheno.selected$pheno)
  {
    data <- query.pheno(phenoname)
    facet <- phenoname
    t <- pheno.selected[pheno.selected$pheno ==facet,]$search.term
    if (!is.na(t)) {
      facet <- paste0(facet," - ", t)
    }
    append <- data.frame(pred_perf_r2 = data$pred_perf_r2,
                         zscore = data$zscore,
                         pval = data$pval,
                         the_facet = facet)
    d <- rbind(d, append)
  }
  return(d)
}

find_snps_for_genes_more_significant_than_gwas <- function(gwas, data, pheno_tag, window_semi=1e6) {
  m_results_p <- data %>% 
    filter(phenotype == pheno_tag) %>% 
    select(phenotype, gene, gene_name, tissue, pval)
  
  m_results_p <- m_results_p %>% inner_join(gencode, by = c("gene" = "gene_key", "gene_name"="gene_name"))
  #explicit logic asdplyr doesnt handel inequality join
  # and sqldf is too solow
  n <- length(m_results_p$gene)
  m <- m_results_p %>% filter(pval < 0.05/n) %>% group_by(gene) %>% top_n(n = 1, wt = -pval)
  m <- m %>%  
    mutate(window_start = ifelse(start_location-window_semi<0, 0, start_location-window_semi),
           window_end = end_location+window_semi)
  top_snp <- data.frame()
  for (i in 1:length(m$gene)) {
    the_gene <- m$gene[i]
    row <- m %>% filter(gene == the_gene)
    g_chr <- row$chr[1]
    w_s <- row$window_start[1]
    w_e <- row$window_end[1]
    snps <- cad_gwas %>% filter(chrx==g_chr, w_s<=start_position, start_position<=w_e) 
    snps$gene <- the_gene
    t <- snps[which.min(snps$pvalue),]
    if (length(t$rsid) <1) next
    if (t$pvalue[1] < m$pval[i]) next
    top_snp <- rbind(top_snp, data.frame(gene=snps$the_gene, rsid=snps$rsid, pvalue=snps$pvalue))
    print(paste0(i," ", m$gene_name[i]))
  }
  m <- m %>% inner_join(top_snp, by="gene") %>%
    select(phenotype, gene, tissue, gene_name, pval, top_rsid, top_snp_pvalue)
  return(m)
}

s <- find_snps_for_genes_more_significant_than_gwas(cad_gwas, d, "CARDIoGRAM_C4D_CAD_ADDITIVE" )
