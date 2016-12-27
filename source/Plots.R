#!/usr/bin/env Rscript

########################################################################################################################
# Handy plotting and helper functions.

library(ggplot2)
library(qqman)

do_manhattan_plot_from_data <- function(data, file_prefix) {
    c <- data[complete.cases(data),]
    if(!("pval" %in% colnames(c)))
    {
        c$pval <- 2*pnorm(-abs(c$zscore))
    }

    nn = length(data$zscore)
    p_b = 0.05/nn

    image <- paste(file_prefix, "-manhattan.png", sep="")
    png(filename=image,width=1024,height=768)
    manhattan(c, chr="chr", bp="base_position", p="pval", snp="gene_name",
            suggestiveline=-log10(p_b), genomewideline=FALSE, annotatePval=p_b, annotateTop=FALSE)
#    manhattan(c, chr="chr", bp="base_position", p="pval", snp="gene", annotatePval=1e-5, annotateTop=FALSE)
    dev.off()
}