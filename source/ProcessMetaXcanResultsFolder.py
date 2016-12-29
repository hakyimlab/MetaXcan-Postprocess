#! /usr/bin/env python
__author__ = 'heroico'

import os
import re
import logging
import numpy
import Gencode
import Logging
import Utilities
import pandas

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
dplyr = importr('dplyr', on_conflict="warn")
#import rpy2.robjects.lib.dplyr as dplyr
robjects.r.source("Plots.R")
robjects.r.source("PlotsQQ.R")

from rpy2.robjects import numpy2ri
numpy2ri.activate()

def load_metaxcan_result(path):
    d = pandas.read_csv(path)
    d = d.iloc[numpy.isfinite(d.zscore.values)]
    d = d[["gene_name", "zscore", "pvalue"]]

    d = robjects.r['data.frame'](
        gene_name = robjects.StrVector(d["gene_name"].values.tolist()),
        zscore = robjects.FloatVector(d["zscore"].values.tolist()),
        pvalue = robjects.FloatVector(d["pvalue"].values.tolist()),
        stringsAsFactors=False
    )

    return d

def annotation(path, only_numeric_chromosomes=True):
    CHROMOSOME = Gencode.GFTF.K_CHROMOSOME
    START_LOCATION = Gencode.GFTF.K_START_LOCATION
    GENE_NAME = Gencode.GFTF.K_GENE_NAME
    gencode = Gencode.load_gene_annotation(path)

    gencode.chromosome = gencode[CHROMOSOME].str.replace("chr", "")
    if only_numeric_chromosomes:
        gencode = gencode[gencode[[CHROMOSOME]].apply(lambda x: x[0].isdigit(), axis=1)]

    gencode = robjects.r['data.frame'](
        chromosome=robjects.StrVector(gencode[CHROMOSOME].values.tolist()),
        start_location=robjects.IntVector(gencode[START_LOCATION].values.tolist()),
        gene_name = robjects.StrVector(gencode[GENE_NAME].values.tolist()),
        stringsAsFactors=False)

    return gencode

def process(gene_annotation, in_path, out_prefix):
    logging.info("Processing %s", in_path)
    d = load_metaxcan_result(in_path)

    qq_path = out_prefix+"_qqplot.png"
    robjects.r.do_qq_unif_plot(d, qq_path, "QQ Plot")

    m = dplyr.inner_join(d, gene_annotation, by="gene_name")
    manhattan_path = out_prefix+"_manhattan.png"
    robjects.r.do_manhattan_plot_from_data(m, manhattan_path)

def run(args):
    logging.info("Loading gencode")
    gene_annotation = annotation(args.gencode_file)

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    pattern = re.compile(args.pattern) if args.pattern else None
    files = Utilities.filesFromFolder(args.input_folder, pattern)
    for f in files:
        in_path = os.path.join(args.input_folder, f)
        out_path = os.path.join(args.output_folder, f)
        out_prefix = os.path.splitext(out_path)[0]
        process(gene_annotation, in_path, out_prefix)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Batch processing of metaxcan results in a folder')

    parser.add_argument("--input_folder",
                        help="name of folder with zscore files",
                        default=None)

    parser.add_argument("--pattern",
                        help="name pattern to select files",
                        default=None)

    parser.add_argument("--output_folder",
                        help="where to output stuff",
                        default="results/images")

    parser.add_argument("--gencode_file",
                        help="path of gencode gene information digest",
                        default="data/gencode.v18.genes.patched_contigs.summary.protein")

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--throw",
                        action="store_true",
                        help="Throw exception on error",
                        default=False)

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    run(args)