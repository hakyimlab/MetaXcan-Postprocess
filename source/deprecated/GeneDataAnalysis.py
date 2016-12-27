#! /usr/bin/env python

import logging
import os

from IPython import embed
import numpy

import metax.WeightDBUtilities as WeightDBUtilities
import metax.Logging as Logging
import metax.Utilities as Utilities
import metax.KeyedDataSet as KeyedDataSet
import metax.MatrixUtilities as MatrixUtilities
import metax.ZScoreCalculation as ZScoreCalculation
import metax.MethodGuessing as MethodGuessing



def get_gene_data(gene_name, weight_db_logic, covariance_contents, betas):
    g = None
    gs = [gene_entry for gene, gene_entry in weight_db_logic.gene_data_for_gene.iteritems() if gene_entry.gene_name == gene_name]
    if len(gs): g = gs[0]
    else: raise RuntimeError("No valid gene data found")

    weights = weight_db_logic.weights_by_gene[g.gene]

    covariance_matrix, valid_rsids = covariance_contents[g.gene]

    beta_sets = None
    for k,b in betas.iteritems():
        for rsid in weights:
            snps = b["rsid"].values_by_key
            if rsid in snps:
                beta_sets = b
    return g, weights, covariance_matrix, valid_rsids, beta_sets


def run(args):
    logging.info("Loading weight db")
    weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(args.weight_db_path)

    logging.info("Loading covariance file")
    covariance_contents = MatrixUtilities.loadMatrixFromFile(args.covariance)

    logging.info("Choosing method")
    beta_contents = Utilities.contentsWithPatternsFromFolder(args.beta_folder, [])
    zscore_calculation, normalization = MethodGuessing.chooseZscoreSchemeFromFiles(args.beta_folder, beta_contents, covariance_contents, weight_db_logic)

    logging.info("Processing")
    betas = {}
    for content in beta_contents:
        logging.info("Loading betas")
        beta_path = os.path.join(args.beta_folder, content)
        beta_sets = KeyedDataSet.KeyedDataSetFileUtilities.loadDataSetsFromCompressedFile(beta_path, header="")
        beta_sets = {set.name: set for set in beta_sets}
        betas[content] = beta_sets

    if args.gene_name:
        try:
            gene_data, weights, covariance_matrix, valid_rsids, beta_sets = get_gene_data(args.gene_name, weight_db_logic, covariance_contents, betas)
            weight_values, variances = ZScoreCalculation.preProcess(covariance_matrix, valid_rsids, weights, beta_sets)
            if args.interactive:
                embed()
            logging.info("Processed gene data")
        except Exception as e:
            logging.info("Couldn't get gene data")
            embed()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Gather MetaXcan related data.')

    parser.add_argument("--beta_folder",
                        help="name of weight db in data folder",
                        default="intermediate/beta/ANGST")

    parser.add_argument("--weight_db_path",
                        help="name of weight db in data folder",
                        default="data/dbs/TW_Brain-Anteriorcingulatecortex-BA24_0.5.db")

    parser.add_argument("--covariance",
                        help="name of folder containing covariance data",
                        default="intermediate/cov/TW_Brain-Anteriorcingulatecortex-BA24.txt.gz")

    parser.add_argument("--gene_name",
                        help="name of gene",
                        default=None)

    parser.add_argument("--output_file",
                        help="name of output file",
                        default="forensics/gene_forensics.csv")

    parser.add_argument("--split_output_prefix",
                        help="name of output file",
                        default="forensics/OR4C46")

    parser.add_argument("--interactive",
                    help="Wether input files are gzip compressed file",
                    action="store_true",
                    default=False)

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    run(args)