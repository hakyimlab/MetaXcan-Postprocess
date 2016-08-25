
find.clinvar.genes = function(clinvar, disname) unique( clinvar$GeneSymbol[grep(disname,clinvar$DiseaseName,ignore.case=T )] )

build_qqunif_data <- function(data, clinvar, pheno.selected, genesetname) {
  d1 <- data.frame()
  d2 <- data.frame()
  
  parse_pheno <- function(in.df, facet){
    the_facet <- facet
    p_val <- 2*pnorm(-abs(in.df$zscore))
    y <- -sort(log10(p_val))
    nn <- length(y)
    x <- -log10((1:nn)/(nn+1))
    
    c95 <- rep(0,nn)
    c05 <- rep(0,nn)
    for(j in 1:nn)
    {
      c95[j] <- qbeta(0.95,j,nn-j+1)
      c05[j] <- qbeta(0.05,j,nn-j+1)
    }
    y95 <- -log10(c95)
    y05 <- -log10(c05)
    b <- -log10(0.05/nn) #bonferroni
    
    y_fdr005 <- x -log10(0.05)
    y_fdr010 <- x - log10(0.1)
    y_fdr025 <- x - log10(0.25)
    
    built <- data.frame(x=x, y=y, y95=y95, y05=y05, b=b, 
                        y_fdr005 = y_fdr005, y_fdr025 = y_fdr025, y_fdr010 = y_fdr010,
                        the_facet=the_facet, stringsAsFactors = FALSE)
    return(built)    
  }
  
  for(phenoname in pheno.selected$pheno)
  {
    disname = pheno.selected$search.term[pheno.selected$pheno == phenoname]
    if (is.na(disname)) {
      next;
    }
    geneset = find.clinvar.genes(clinvar, disname)
    #geneset = clinvargene$GeneSymbol
    
    d <- data[data$phenotype == phenoname,]
    ind <- d$gene_name %in% geneset
    ngenes <- length(unique(d$gene_name[ind]))
    phenotype <- unique(d$phenotype)
    phenotype.short <- pheno.selected$pheno.short[pheno.selected$pheno==phenoname]
    title = paste0(phenotype.short,' ', ', all genes vs (', length(geneset),') ', disname,' genes')
    
    built1 <- parse_pheno(d, title)
    
    #repeat only for selected genes
    d <- d[ind,]
    built2 <- parse_pheno(d, title)
    
    if(nrow(built2)>1)
    {
      d1 <- rbind(d1, built1)
      d2 <- rbind(d2, built2)
    }
  }
  return(list(d1=d1,d2=d2))
}

build_decoration_data <- function(data) {
  r <- data.frame()
  x <- d$x
  y_fdr005 <- x -log10(0.05)
  r <- rbind(r, data.frame(x=x, y=y_fdr005, tag=rep('FDR = 0.05', length(x))), the_facet = d$the_facet)
  y_fdr010 <- x - log10(0.1)
  r <- rbind(r, data.frame(x=x, y=y_fdr010, tag=rep('FDR = 0.10', length(x))), the_facet = d$the_facet)
  y_fdr25 <- x - log10(0.25)
  r <- rbind(r, data.frame(x=x, y=y_fdr25, tag=rep('FDR = 0.25', length(x))), the_facet = d$the_facet)
  return(r)
}

qq_unif_faceted_comparison_plot <- function(d1, d2, title=NULL, columns=2) {
  p <- ggplot(d1) +
    theme_bw() +
    theme(plot.title = element_text(lineheight=3, face="bold", size=25)) +
    xlab(expression(Expected~~-log[10](italic(p)))) +
    ylab(expression(Observed~~-log[10](italic(p)))) +
    geom_point(data = d1, aes(x=x, y=y), colour = 'black', shape=1) +
    geom_point(data = d2, aes(x=x, y=y), colour = 'blue', shape=1) +
    facet_wrap(~the_facet, scales="free",ncol=columns)
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  p <- p + geom_abline(intercept=0, slope=1, colour='gray')

  max_x <- max(d1$x)
  x <- seq(0, max_x)
  y005 <- seq(-log10(0.05), max_x-log10(0.05))
  
  decoration <- build_decoration_data(d1)
  p <- p +
    geom_line(mapping = aes(x=x, y=y_fdr005, colour = 'FDR = 0.05'), show_guide = TRUE) +
    geom_line(mapping = aes(x=x, y=y_fdr010, colour = 'FDR = 0.10'), show_guide = TRUE) +
    geom_line(mapping = aes(x=x, y=y_fdr025, colour = 'FDR = 0.25'), show_guide = TRUE) +
    scale_colour_manual(labels = c('FDR = 0.05', 'FDR = 0.10', 'FDR = 0.25'),
                        breaks = c('FDR = 0.05', 'FDR = 0.10', 'FDR = 0.25'),
                        values = c("red", "orange", "yellow")) +
    theme(legend.position="bottom") + 
    theme(legend.text=element_text(size=25)) +
    theme(legend.title=element_text(size=25)) +
    theme(legend.key.size = unit(50, "points")) +
    guides(color=guide_legend("FDR levels"))
  
    
  p <- p + geom_abline(data=d1,aes(intercept=b, slope=0), colour='black') ## bonferroni
  #
  p <- p + geom_line(data=d1, aes(x=x, y=y95), colour='gray')
  p <- p + geom_line(data=d1, aes(x=x, y=y05), colour='gray')
  return(p)
}


