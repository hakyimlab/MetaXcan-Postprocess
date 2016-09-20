#! /usr/bin/env Rscript
library(argparse)
library(dplyr)
library(ggplot2)

parser <- ArgumentParser(description="Compare matching metaxcan results folders")
parser$add_argument('--stats_file',
                    help='path to stats',
                    default='data/paper_stats.csv')

parser$add_argument('--beta_lambda_file',
                    help='path to lambda',
                    default='results/beta_lambda.csv')

parser$add_argument('--merged_file',
                    help='results/lambda_merged.csv',
                    default='results/merged.csv')

parser$add_argument('--plot',
                    help='image output',
                    default='results/lambda.png')

arguments <- parser$parse_args(commandArgs(TRUE))

l <- read.csv(arguments$stats_file, stringsAsFactors = FALSE)
l[l$InflationLambda == 0,]$InflationLambda <- NA
r <- read.csv(arguments$beta_lambda_file, stringsAsFactors = FALSE)
r <- rename(r, Tag = pheno)
r <- rename(r, BetaLambda = lambda)
j <- left_join(l, r, by = "Tag")

p <- ggplot( j, aes(x=InflationLambda, y=BetaLambda)) +
        theme_bw() +
        theme(plot.title = element_text(lineheight=3, face="bold")) +
        xlab("MetaXcan Inflation") +
        ylab("GWAS Inflation") +
        geom_point(colour = 'black', shape=1) +
        geom_abline(intercept=0, slope=1, colour="grey49")

png(filename=arguments$plot,width=1024,height=1024)
print(p)
dev.off()

write.csv(j, arguments$merged_file)


