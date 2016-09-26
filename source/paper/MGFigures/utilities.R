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

split_into_lines <- function(strings, size) {
  u <- unique(strings)
  l <- list()
  ppaste <- function(a, b, p) {
    if (nchar(a) > 0) {
      a <- paste0(a,"_")
    }
    a <- paste0(a,b,p)
    a
  }
  for (word in u) {
    comps <- strsplit(word, "_")[[1]]
    line <- ""
    buff <- ""
    for (comp in comps) {
      if (nchar(buff) + nchar(comp) > 19 && length(comps) > 1) {
        line <- ppaste(line, buff, "\n")
        buff <- comp
      } else {
        buff <- ppaste(buff, comp, "")
      }
    }
    line <- ppaste(line, buff, "")
    l[[word]] <- line
  }
  return(l)
}
