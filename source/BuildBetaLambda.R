#! /usr/bin/env Rscript
library(argparse)
library(dplyr)

#Math
calc_lambda <- function(p.vec){
  #http://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html
  chisq <- qchisq(1-na.omit(p.vec),1)
  lambda <- median(chisq)/qchisq(0.5,1)
  return(lambda)
}

gc_adjust <- function(z.vec,p.vec){
  lambda <- calc_lambda(p.vec)
  #print("Inflation factor: " %&% lambda)
  z.adj.vec <- z.vec/lambda #z.vec/sqrt(lambda)  #
  chi.vec <- z.adj.vec^2
  p.adj.vec <- pchisq(chi.vec,1,lower.tail=FALSE) #pnorm(-abs(z.adj.vec))*2
  return(list(z.adj.vec,p.adj.vec))
}

#file system bureaucracy
get_path_df <- function(folder, pheno_remove_prefix) {
    folder_names <- list.files(folder)

    parse_names <- function(names, split, tissue_pre_addition, pheno_pre_removal) {
        comps <- strsplit(names, split)
        parse <- function(x, i, addition) {
            if (length(x) <= 1) {
                return("")
            } else {
                return(sprintf("%s%s", c(addition), x[[i]]))
            }
        }
        phenos <- unlist(lapply(comps, function(x) parse(x,1, "")))
        phenos <- strsplit(phenos, pheno_pre_removal)
        phenos <- unlist(lapply(phenos, function(x) parse(x, 2, "")))

        tissues <- unlist(lapply(comps, function(x) parse(x,2, tissue_pre_addition)))

        return( list(pheno=phenos, tissue=tissues) )
    }

    comps_tw <- parse_names(folder_names, "_TW", "TW", pheno_remove_prefix)
    comps_dgn <- parse_names(folder_names, "_DGN", "DGN", pheno_remove_prefix)
    comps_cross <- parse_names(folder_names, "_Cross", "Cross",pheno_remove_prefix)
    merge <- function(a, b, c) {
        r <- rep(a)
        r <- pmax(r, b)
        r <- pmax(r, c)
        return(r)
    }
    pheno <- merge(comps_tw$pheno, comps_dgn$pheno, comps_cross$pheno)
    tissue <- merge(comps_tw$tissue, comps_dgn$tissue, comps_cross$tissue)

    path <- sprintf("%s/%s", folder, folder_names)

    d <- data.frame(path=path, pheno = pheno, tissue = tissue, stringsAsFactors = FALSE)
    return(d)
}

get_pvalues <- function(file) {
    df <- read.delim(file, sep = " ")
    pvalue <- 2*pnorm(-abs(df$beta_z))
    return(pvalue)
}

get_lambdas <- function(df) {
    phenos <- distinct(df, pheno)$pheno
    lambdas <- numeric()
    for (the_pheno in phenos) {
        print(paste0("Processing pheno: ",the_pheno,"..."))
        filtered <- df %>% filter(pheno == the_pheno)
        p <- numeric()
        for (folder in filtered$path) {
            files <- list.files(folder)
            for (file in files) {
                path <- sprintf("%s/%s",folder,file)
                print(paste0("Processing ", path, " ..."))
                p <- rbind(p, get_pvalues(path))
            }
        }
        lambdas <- rbind(lambdas, calc_lambda(p))
    }
    result <- data.frame(pheno = phenos, lambda = lambdas)
    return(result)
}

parser <- ArgumentParser(description="Compare matching metaxcan results folders")
parser$add_argument('--folder',
                    help='path to old results',
                    default='data/betas')

parser$add_argument('--beta_prefix',
                    help='folder prefix',
                    default='beta_')

parser$add_argument('--result',
                    help='result',
                    default='results/beta_lambda.csv')


arguments <- parser$parse_args(commandArgs(TRUE))

d <- get_path_df(arguments$folder, arguments$beta_prefix)

l <- get_lambdas(d)

write.csv(l, arguments$result)

