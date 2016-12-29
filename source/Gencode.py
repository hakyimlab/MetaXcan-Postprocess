__author__ = 'heroico'
#trimmed from PredictDBAnalysis/gencode_input

import csv
import gzip
import pandas

K_NOT_GENES = ["transcript","exon","CDS","UTR","start_codon","stop_codon","Selenocysteine"];

# look at gencode http://www.gencodegenes.org/data_format.html
class GFTF:
    """gencode file table format"""
    CHROMOSOME = 0
    ANNOTATION_SOURCE = 1
    FEATURE_TYPE = 2
    N_START_LOCATION = 3
    N_END_LOCATION = 4
    SCORE = 5
    GENOMIC_STRAND = 6
    GENOMIC_PHASE = 7
    KEY_VALUE_PAIRS = 8

    #there are several other key-value pairs but we are concerned with these
    K_CHROMOSOME = "chromosome"
    K_START_LOCATION = "start_location"
    K_GENE_ID = "gene_id"
    K_TRANSCRIPT_ID = "transcript_id"
    K_GENE_TYPE = "gene_type"
    K_GENE_STATUS = "gene_status"
    K_GENE_NAME = "gene_name"
    K_TRANSCRIPT_TYPE = "transcript_type"
    K_TRANSCRIPT_STATUS = "transcript_status"
    K_TRANSCRIPT_NAME = "transcript_name"
    K_EXON_NUMBER = "exon_number"
    K_EXON_ID = "exon_id"
    K_LEVEL = "level"

    #some are missing
    K_TAG = "tag"

def key_value_pairs_to_dict(key_value_pairs):
    d = {}
    key, value = None, None
    for i, string in enumerate(key_value_pairs):
        if key is None:
            key = string
        elif value is None:
            value = string.translate(None, '\'";\n')
            d[key] = value
            key = None
            value = None
    return d

def load_gene_annotation(path, only_genes=True):
    gene_id = []
    gene_name = []
    chromosome =[]
    start_location = []

    with gzip.open(path) as file:
        for i,line in enumerate(file):
            comps = line.strip().split()

            if "##" in comps[0]:
                continue

            if only_genes and comps[GFTF.FEATURE_TYPE] in K_NOT_GENES:
                continue

            key_value_pairs = [x.translate(None, "'\";") for x in comps[GFTF.KEY_VALUE_PAIRS:]]
            key_value_pairs = key_value_pairs_to_dict(key_value_pairs)

            chromosome.append(comps[GFTF.CHROMOSOME])
            start_location.append(comps[GFTF.N_START_LOCATION])
            gene_id.append(key_value_pairs[GFTF.K_GENE_ID])
            gene_name.append(key_value_pairs[GFTF.K_GENE_NAME])


    r = {GFTF.K_CHROMOSOME:chromosome,
         GFTF.K_START_LOCATION:start_location,
         GFTF.K_GENE_ID:gene_id,
         GFTF.K_GENE_NAME:gene_name
        }

    r = pandas.DataFrame(r)
    r = r.drop_duplicates()
    return r