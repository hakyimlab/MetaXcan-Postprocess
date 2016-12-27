#! /usr/bin/env python
__author__ = 'heroico'

import os
import re
import logging
import rpy2.robjects as robjects
from subprocess import  call
import Logging

def run(args):
    pass

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Batch processing of metaxcan results in a folder')

    parser.add_argument("--input_folder",
                        help="name of folder with zscore files",
                        default="results_ew/results")

    parser.add_argument("--pattern",
                        help="name pattern to select files",
                        default="AMD.*")

    parser.add_argument("--output_folder",
                        help="where to output stuff",
                        default="results_ew/images_AMD")

    parser.add_argument("--gene_digest_file",
                        help="path of gene information digest",
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