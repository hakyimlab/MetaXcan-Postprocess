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