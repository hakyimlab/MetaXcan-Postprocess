#!/usr/bin/env Rscript

########################################################################################################################
# Handy plotting and helper functions.

library(ggplot2)
library(qqman)

do_manhattan_plot_from_data <- function(data, path) {
    c <- data[complete.cases(data),]
    if(!("pvalue" %in% colnames(c)))
    {
        c$pvalue <- 2*pnorm(-abs(c$zscore))
    }

    c$pvalue <- pmax(c$pvalue, 1e-30)

    nn = length(data$zscore)
    p_b = 0.05/nn

    c$chromosome <- as.numeric(c$chromosome)

    png(filename=path,width=1024,height=768)
    manhattan(c, chr="chromosome", bp="start_location", p="pvalue", snp="gene_name",
            suggestiveline=-log10(p_b), genomewideline=FALSE, annotatePval=p_b, annotateTop=FALSE)
#    manhattan(c, chr="chr", bp="base_position", p="pval", snp="gene", annotatePval=1e-5, annotateTop=FALSE)
    dev.off()
}