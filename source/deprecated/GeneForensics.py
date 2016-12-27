#! /usr/bin/env python
__author__ = 'heroico'

########################################################################################################################
# Gathers some figures that go into MEtaXcan's calculation for a given gene

import logging
import os
import numpy

import metax.WeightDBUtilities as WeightDBUtilities
import metax.Logging as Logging
import metax.Utilities as Utilities
import metax.KeyedDataSet as KeyedDataSet
import metax.MatrixUtilities as MatrixUtilities
import metax.ZScoreCalculation as ZScoreCalculation
import metax.MethodGuessing as MethodGuessing

def get_gene_data(weight_db_logic, gene_name):
    g = None
    for gene, gene_entry in weight_db_logic.gene_data_for_gene.iteritems():
        if gene_entry.gene_name == args.gene_name:
            g = gene_entry
            break
    if g == None:
        raise RuntimeError("No valid gene data found")
    return g

def build_split_output(args, covariance_matrix, valid_rsids, weights, zscore_calculation, beta_sets, weight_values, variances):
    logging.info("Building split output")
    with open(args.split_output_prefix+"_snps.txt", "w") as file:
        file.write("rsid, w_l, z_l, var_l\n")
        for i in xrange(0, len(valid_rsids)):
            rsid = valid_rsids[i]
            w = weights[rsid] if rsid in weights else 0
            b = beta_sets["beta"].values_by_key[rsid] if rsid in beta_sets["beta"].values_by_key else 0
            b_z = beta_sets["beta_z"].values_by_key[rsid] if rsid in beta_sets["beta_z"].values_by_key else 0
            v = variances[rsid] if rsid in variances else 0
            file.write(",".join([rsid, str(w.weight), str(b_z), str(v)])+"\n")

    with open(args.split_output_prefix+"_cov.txt", "w") as file:
        file.write(",".join(valid_rsids)+"\n")
        for i in xrange(0, len(valid_rsids)):
            row = map(str,list(covariance_matrix[i]))
            file.write(",".join(row)+"\n")


def build_output(args, covariance_matrix, valid_rsids, weights, zscore_calculation, beta_sets, weight_values, variances):
    logging.info("Building forensics")
    folder = os.path.split(args.output_file)[0]
    if len(folder) and not os.path.exists(folder):
        os.makedirs(folder)
    with open(args.output_file, "w") as file:
        header = [""]
        header = header + valid_rsids + ["","w","b", "b_z", "v", "\n"]
        file.write(",".join(header))
        n = len(header)
        for i in xrange(0, len(valid_rsids)):
            rsid = valid_rsids[i]
            w = weights[rsid] if rsid in weights else 0
            b = beta_sets["beta"].values_by_key[rsid] if rsid in beta_sets["beta"].values_by_key else 0
            b_z = beta_sets["beta_z"].values_by_key[rsid] if rsid in beta_sets["beta_z"].values_by_key else 0
            v = variances[rsid] if rsid in variances else 0
            row = [rsid] + map(str,list(covariance_matrix[i])) + ["", str(w.weight), str(b), str(b_z), str(v)]
            file.write(",".join(row)+"\n")

        blank_row = ",".join(["" for i in xrange(0,n)])+"\n"
        file.write(blank_row)
        file.write(blank_row)
        file.write(blank_row)

        eigenvalues, eigenvectors = numpy.linalg.eigh(covariance_matrix)
        eigen_row = ["eigenvalues"] + map(str, list(eigenvalues)) + ["" for i in xrange(0, n-len(eigenvalues)-2)]
        file.write(",".join(eigen_row)+"\n")

        dot_product = numpy.dot(numpy.dot(numpy.transpose(weight_values), covariance_matrix), weight_values)
        dot_row = ["dot_product", str(dot_product)] + ["" for i in xrange(0,n-2)]
        file.write(",".join(dot_row)+"\n")

        det = numpy.linalg.det(covariance_matrix)
        det_row = ["det", str(det)] + ["" for i in xrange(0,n-2)]
        file.write(",".join(det_row)+"\n")

        eigenmax = numpy.amax(eigenvalues)
        eigenmax_row = ["max eigenvalue", str(eigenmax)] + ["" for i in xrange(0,n-2)]
        file.write(",".join(eigenmax_row)+"\n")

        eigenmin = numpy.amin(eigenvalues)
        eigenmin_row = ["min eigenvalue", str(eigenmin)] + ["" for i in xrange(0,n-2)]
        file.write(",".join(eigenmin_row)+"\n")

        pre_zscore, n_snps, VAR_g, effect_size = zscore_calculation(args.gene_name, weights, beta_sets, covariance_matrix, valid_rsids)
        zscore_row = ["zscore", str(pre_zscore)]+ ["" for i in xrange(0,n-2)]
        file.write(",".join(zscore_row)+"\n")
        VAR_g_row = ["VAR_g", str(VAR_g)]+ ["" for i in xrange(0,n-2)]
        file.write(",".join(VAR_g_row)+"\n")
        n_snps_row = ["n_snps", str(n_snps)]+ ["" for i in xrange(0,n-2)]
        file.write(",".join(n_snps_row)+"\n")

def run(args):
    logging.info("Loading weight db")
    weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(args.weight_db_path)
    gene_data = get_gene_data(weight_db_logic, args.gene_name)
    weights = weight_db_logic.weights_by_gene[gene_data.gene]

    logging.info("Loading covariance file")
    covariance_contents = MatrixUtilities.loadMatrixFromFile(args.covariance)
    covariance_matrix, valid_rsids = covariance_contents[gene_data.gene]

    logging.info("Choosing method")
    beta_contents = Utilities.contentsWithPatternsFromFolder(args.beta_folder, [])
    zscore_calculation, normalization = MethodGuessing.chooseZscoreSchemeFromFiles(args.beta_folder, beta_contents, covariance_contents, weight_db_logic)

    logging.info("Processing")
    for content in beta_contents:
        logging.info("Loading betas")
        beta_path = os.path.join(args.beta_folder, content)
        beta_sets = KeyedDataSet.KeyedDataSetFileUtilities.loadDataSetsFromCompressedFile(beta_path, header="")
        beta_sets = {set.name:set for set in beta_sets }
        weight_values, variances = ZScoreCalculation.preProcess(covariance_matrix, valid_rsids, weights, beta_sets)

        for rsid in valid_rsids:
            if rsid in beta_sets["beta_z"].values_by_key:
                build_output(args, covariance_matrix, valid_rsids, weights, zscore_calculation, beta_sets, weight_values, variances)
                build_split_output(args, covariance_matrix, valid_rsids, weights, zscore_calculation, beta_sets, weight_values, variances)
                break


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build betas from GWAS data.')

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
                        default="OR4C46")

    parser.add_argument("--output_file",
                        help="name of output file",
                        default="forensics/gene_forensics.csv")

    parser.add_argument("--split_output_prefix",
                        help="name of output file",
                        default="forensics/OR4C46")

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    run(args)