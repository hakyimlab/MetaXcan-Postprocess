#!/usr/bin/env Rscript

########################################################################################################################
# Handy plotting and helper functions for producing QQ Plots

library(dplyr)
library(tidyr)

build_multi_qqunif_data <- function(data) {
  parse_pheno <- function(in.df, facet){
    the_facet <- facet
    if("pvalue" %in% colnames(in.df)) {
      p_val <- in.df$pvalue
    } else {
      p_val <- 2*pnorm(-abs(in.df$zscore))
    }

    y <- -sort(log10(p_val))
    y <- pmin(y, 30) #upper threshold value
    nn <- length(y)
    x <- -log10((1:nn)/(nn+1))
    b <- -log10(0.05/nn) #bonferroni

    y_fdr005 <- x -log10(0.05)
    y_fdr010 <- x - log10(0.1)
    y_fdr025 <- x - log10(0.25)

    built <- data.frame(x=x, y=y, b=b,
                        y_fdr005 = y_fdr005, y_fdr025 = y_fdr025, y_fdr010 = y_fdr010,
                        the_facet=the_facet, stringsAsFactors = FALSE)
    return(built)
  }

  tissues <- as.character(unique(data$tissue))
  phenos <- as.character(unique(data$pheno))
  r <- data.frame()
  for(t in tissues)
  {
    td <- data %>% filter(tissue == t)
    for (p in phenos) {
      ptd <- td %>% filter(pheno == p)
      built <- parse_pheno(ptd, t)
      built$pheno <- p
      r <- rbind(r, built)
    }
  }
  return(r)
}

qq_unif_faceted_comparison_plot <- function(data, title=NULL, columns=2) {
  p <- ggplot() +
    theme_bw() +
    theme(plot.title = element_text(lineheight=3, face="bold", size=25)) +
    xlab(expression(Expected~~-log[10](italic(p)))) +
    ylab(expression(Observed~~-log[10](italic(p)))) +
    facet_wrap(~the_facet, scales="free",ncol=columns)

  values <- c('FDR = 0.05' = "red", 'FDR = 0.10' = "orange", 'FDR = 0.25' = "yellow")
  palette <- c("blueviolet", "cyan4", "chartreuse4", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  n <- unique(data$the_facet)
  for (i in 1:length(n)) {
    name <- n[i]
    values[name] <- palette[i]
  }

  p <- p +
    scale_colour_manual(values = values) +
    theme(legend.position="bottom") +
    theme(legend.text=element_text(size=20)) +
    theme(legend.title=element_text(size=20)) +
    theme(legend.key.size = unit(40, "points")) +
    theme(strip.text.x = element_text(size=18, face="bold"),
          strip.background = element_rect(fill="white")) +
    guides(color=guide_legend("Colors :"))

  #the_facet <- data$the_facet[1]
  #d <- data %>% filter(the_facet== the_facet)
  p <- p + geom_line(data = data, mapping = aes(x=x, y=y_fdr005, colour = 'FDR = 0.05'), show.legend = TRUE) +
    geom_line(data = data, mapping = aes(x=x, y=y_fdr010, colour = 'FDR = 0.10'), show.legend = TRUE) +
    geom_line(data = data, mapping = aes(x=x, y=y_fdr025, colour = 'FDR = 0.25'), show.legend = TRUE) +
    geom_point(data = data, mapping=aes(x=x, y=y, group=the_facet, colour=the_facet), size=2)

  p <- p + geom_abline(data=d,mapping= aes(intercept=b, slope=0), colour='black')
  p <- p + geom_abline(intercept=0, slope=1, colour='gray')

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

build_qqunif_data <- function(data, label=NULL) {
  if("pvalue" %in% colnames(data)) {
    p_val <- data$pvalue
  } else {
    p_val <- 2*pnorm(-abs(data$zscore))
  }

  y <- -sort(log10(p_val))
  y <- pmin(y, 30) #upper threshold value
  nn <- length(y)
  x <- -log10((1:nn)/(nn+1))
  b <- -log10(0.05/nn) #bonferroni

  y_fdr005 <- x -log10(0.05)
  y_fdr010 <- x - log10(0.1)
  y_fdr025 <- x - log10(0.25)

  built <- data.frame(x=x, y=y, b=b,
                      y_fdr005 = y_fdr005, y_fdr025 = y_fdr025, y_fdr010 = y_fdr010,
                     stringsAsFactors = FALSE)
  if (!is.null(label)) {
    built$label <- label;
  }
  return(built)
}

build_qqunif_data_using_label <-function(data) {
  results <- data.frame()
  labels <- unique(data$label)
  for (l in labels) {
    d <- data %>% filter(label == l)
    q <- build_qqunif_data(d, l)
    results <- rbind(results, q)
  }
  return(results)
}

build_qqunif_data_with_id_label <- function(data) {
  if(!("pvalue" %in% colnames(data))) {
    data$pvalue <- 2*pnorm(-abs(data$zscore))
  }
  sets <- unique(data$label)
  k <- spread(data, label, pvalue)
  first <- sets[1]
  ref <- k[[first]]
  k <- k[order(ref),]

  y <- -log10(ref)
  y <- pmin(y, 30) #upper threshold value
  nn <- length(y)
  k$y <- y
  k$x <- -log10((1:nn)/(nn+1))
  k$b <- -log10(0.05/nn) #bonferroni

  k$y_fdr005 <- x -log10(0.05)
  k$y_fdr010 <- x - log10(0.1)
  k$y_fdr025 <- x - log10(0.25)
  return(k)
}

qq_unif_plot <- function(data, title=NULL, columns=NULL, force_gray_title='GTEx') {
  p <- ggplot() +
    theme_bw() +
    theme(plot.title = element_text(lineheight=3, face="bold", size=15)) +
    xlab(expression(Expected~~-log[10](italic(p)))) +
    ylab(expression(Observed~~-log[10](italic(p))))

  if("facet_w" %in% colnames(data) && !is.null(columns)) {
    p <- p + facet_wrap(~facet_w, scales="free",ncol=columns)
  }

  values <- c('FDR = 0.05' = "red", 'FDR = 0.10' = "orange", 'FDR = 0.25' = "yellow", 'force_gray'='gray')
  palette <- c("blueviolet", "cyan4", "dodgerblue4", "chartreuse4", "#E69F00", "#F0E442", "#999999", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
  if("label" %in% colnames(data)) {
    n <- unique(data$label)
    for (i in 1:length(n)) {
      name <- n[i]
      values[name] <- palette[i]
    }
  }
  p <- p +
    scale_colour_manual(values = values) +
    theme(legend.position="bottom") +
    theme(legend.text=element_text(size=15)) +
    theme(legend.title=element_text(size=12)) +
    theme(legend.key.size = unit(20, "points")) +
    theme(strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_rect(fill="white")) +
    guides(color=guide_legend("Colors :"))

  # if multiple data, show only one as guidelines
  d <- NULL
  if("label" %in% colnames(data)) {
    the_label <- data$label[1]
    d <- data %>% filter(label == the_label)
  } else {
    d <- data
  }
  p <- p + geom_line(data = d, mapping = aes(x=x, y=y_fdr005, colour = 'FDR = 0.05'), show.legend = TRUE) +
    geom_line(data = d, mapping = aes(x=x, y=y_fdr010, colour = 'FDR = 0.10'), show.legend = TRUE) +
    geom_line(data = d, mapping = aes(x=x, y=y_fdr025, colour = 'FDR = 0.25'), show.legend = TRUE)
  if("label" %in% colnames(data)) {
    if ('force_gray' %in% colnames(data)) {
      data$col <- data$label
      data$col[data$force_gray == TRUE] <- 'force_gray'
      d_p <- data %>% filter(data$force_gray != TRUE)
      d_l <- data %>% filter(data$force_gray == TRUE)
      p <- p + geom_line(data = d_l, mapping=aes(x=x, y=y, group=label, colour=factor(col)))
      p <- p + geom_point(data = d_p, mapping=aes(x=x, y=y, group=label, colour=factor(col)), alpha=0.1,size=1)
    } else {
      p <- p + geom_point(data = data, mapping=aes(x=x, y=y, group=label), alpha=0.2,size=1)
    }
  } else {
    p <- p + geom_point(data = data, mapping=aes(x=x, y=y), size=1)
  }

  p <- p + geom_abline(data=d,mapping= aes(intercept=b, slope=0), colour='black')
  p <- p + geom_abline(intercept=0, slope=1, colour='gray')

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

