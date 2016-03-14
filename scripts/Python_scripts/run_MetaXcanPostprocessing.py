# A simple "main()" that loads MetaXcanPostprocessing file 

import sys
from helpers import *
from MetaXcanPostprocessing import * 

if len(sys.argv) != 3:
    hints1 = "Syntax: "
    hints2 = "    python run_MetaXcanPostprocessing.py <gene_list_file> <log_file>"
    print(hints1)
    print(hints2) 
    add_log(hints1)  
    add_log(hints2) 
    exit()

# Start up the logfile
open_log(sys.argv[2])

# Run metaxcan postprocessing 

############################################
#### Part one: Top gene lists with SNPs ####
############################################

# Init class and read .csv file 
snps = TopGeneListSnps(sys.argv[1]) 

# Fetch top gene list from sqlite databases 
snps.fetchTopGeneList(sys.argv[1]) 

# Output data
snps.outputTopGeneList(sys.argv[1]) 

 # Log data 
snps.logTopGeneList(sys.argv[1])





# Close the logfile
finish_log()
