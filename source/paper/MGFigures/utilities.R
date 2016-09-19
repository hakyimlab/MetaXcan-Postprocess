library(stringi)
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

build_allele_key <- function(data.frame) {
  allele_key <- sprintf("%s%s",data.frame$ref_allele, data.frame$eff_allele)
  striHelper <- function(x) stri_c(x[stri_order(x)], collapse = "")
  allele_key <- vapply(stri_split_boundaries(allele_key, type = "character"), striHelper, "")
  return(allele_key)
}
