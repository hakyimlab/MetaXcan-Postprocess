library(dplyr)
source("utilities.R")
#
load_cardiogram_c4d <- function(path) {
  d <- read.delim(path)
  d <- d %>% select(markername,chr,bp_hg19, p_dgc, effect_allele, noneffect_allele)
  d <- d %>% 
    transmute(rsid=markername,
              chr=chr,
              start_position=bp_hg19,
              eff_allele=effect_allele,
              ref_allele=noneffect_allele,
              pvalue=p_dgc)
  d$chrx <- sprintf("chr%d", d$chr)
  #d$allele_key <- build_allele_key(d)
  return(d)
}


load_pgc_bip <- function(path) {
  data <- read.delim(path)
  data <- data %>% 
          select(
            rsid=snpid,
            chr=hg18chr,
            start_position=bp, 
            eff_allele=a1,
            noneffect_allele=a2,
            pvalue=pval)
  data$chrx <- sprintf("chr%d", data$chr)
  return(data)
}
