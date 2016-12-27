__author__ = 'heroico'
#trimmed from PredictDBAnalysis/gencode_input

import csv
import gzip

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
    GENE_ID = "gene_id"
    TRANSCRIPT_ID = "transcript_id"
    GENE_TYPE = "gene_type"
    GENE_STATUS = "gene_status"
    GENE_NAME = "gene_name"
    TRANSCRIPT_TYPE = "transcript_type"
    TRANSCRIPT_STATUS = "transcript_status"
    TRANSCRIPT_NAME = "transcript_name"
    EXON_NUMBER = "exon_number"
    EXON_ID = "exon_id"
    LEVEL = "level"

    #some are missing
    TAG = "tag"