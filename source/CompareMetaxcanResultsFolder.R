#! /usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser(description="Compare matching metaxcan results folders")
parser$add_argument('--first_folder',
                    help='path to old results',
                    default='data/previous')

parser$add_argument('--first_label',
                    help='first label',
                    default='old')

parser$add_argument('--first_folder_selection_pattern',
                    help='String to filter and split files',
                    default='_TW_WholeBlood')

parser$add_argument('--second_folder',
                    help='path to new results',
                    default='data/next')

parser$add_argument('--second_label',
                    help='second label',
                    default='v6p')

parser$add_argument('--second_folder_selection_pattern',
                    help='String to filter and split files',
                    default='_TW_Whole_Blood')

parser$add_argument('--plot_title',
                    help='Cosmetic label',
                    default="Whole blood comparison")

parser$add_argument('--output_path',
                    help='File were to save the results',
                    default='results/comparison.png')

arguments <- parser$parse_args(commandArgs(TRUE))

select_names <- function(folder, pattern) {
    names <- list.files(folder)
    names <- names[grepl(pattern,names)]
    return(names)
}

get_tags <- function(files, pattern) {
    comps <- strsplit(files,pattern)
    first_comps <- lapply(comps, function(x) x[1])
    first_comps <- unlist(first_comps)
    return(first_comps)
}

build_files_in_folder <- function(folder, pattern) {
    names <- select_names(folder, pattern)
    tags <- get_tags(names, pattern)
    f <- paste(strsplit(folder,"/")[[1]],collapse="/")
    files <- unlist(lapply(names, function(x) paste0(f,"/",x) ))
    x <- data.frame(file = files, tag = tags, stringsAsFactors=FALSE)
    return(x)
}

build_files_set <- function(first_folder, first_pattern, second_folder, second_pattern) {
    x <- build_files_in_folder(first_folder, first_pattern)
    y <- build_files_in_folder(second_folder, second_pattern)
    m <- inner_join(x, y, by="tag")
    m <- m[c("file.x", "file.y", "tag")]

    return(m)
}

build_plot_data <- function(file_df) {
    df <- data.frame(x = numeric(), y = numeric(), the_facet = character(), stringsAsFactors=FALSE)
    process_file_row <- function(row){
        x_df <- read.csv(row[1], stringsAsFactors=FALSE)
        y_df <- read.csv(row[2], stringsAsFactors=FALSE)
        m_df <- inner_join(x_df, y_df, by="gene_name")
        df <<- rbind(df, data.frame(x = m_df$zscore.x, y = m_df$zscore.y, the_facet=row[3], stringsAsFactors=FALSE))
    }
    apply(file_df, 1, process_file_row)
    return(df)
}

build_comparison_plot <- function(df, first_label, second_label, title, columns=8) {
    p <- ggplot( df, aes(x=x, y=y)) +
        theme_bw() +
        theme(plot.title = element_text(lineheight=3, face="bold")) +
        xlab(first_label) +
        ylab(second_label) +
        geom_point(colour = 'black', size=0.01) +
        facet_wrap(~the_facet, scales="fixed",ncol=columns) +
        geom_abline(intercept=0, slope=1, colour="grey49")

    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }

    return(p)
}


FIRST_FOLDER <- arguments$first_folder
FIRST_LABEL <- arguments$first_label
FIRST_PATTERN <- arguments$first_folder_selection_pattern

SECOND_FOLDER <- arguments$second_folder
SECOND_LABEL <- arguments$second_label
SECOND_PATTERN <- arguments$second_folder_selection_pattern

TITLE = arguments$TITLE
OUTPUT = arguments$output_path
#Actual script work

file_df <- build_files_set(FIRST_FOLDER, FIRST_PATTERN, SECOND_FOLDER, SECOND_PATTERN)
print(paste0("building plot data for ", length(file_df$file.x), " phenos"))
data <- build_plot_data(file_df)
print("Building plot")
p <- build_comparison_plot(data, FIRST_LABEL, SECOND_LABEL, TITLE)
print("Saving plot")
png(filename=OUTPUT,width=1240,height=2480)
print(p)
dev.off()
print("Done")
