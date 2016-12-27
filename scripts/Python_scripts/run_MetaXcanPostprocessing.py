'''
A simple "main()" that loads functions from MetaXcanPostprocessing.py 
'''

import sys, os, glob, subprocess
from helpers import *
from MetaXcanPostprocessing import * 

if len(sys.argv) != 2:
    hints1 = "Syntax: "
    hints2 = "    python run_MetaXcanPostprocessing.py <project name>"
    exit()

current_path = os.getcwd()

'''
########################
### Open a log fle  ####
########################
'''
# Start up the logfile
open_log(sys.argv[1])


'''
#############################
### Part one: Annotation ####
#############################
'''
# start to log 
pre_message("Part one: Annotating genes ")

os.chdir(current_path)
annotate_metaxcan_result(sys.argv[1])


'''
###########################
#### Part two: QQ-Plot ####
###########################
'''
# start to log 
pre_message("Part two: QQ-Plot ")

os.chdir(current_path)
qqplot(sys.argv[1]) 


'''
####################################
#### Part three: Manhattan-Plot ####
####################################
'''
# start to log 
pre_message("Part three: Manhattan-Plot  ")

os.chdir(current_path)
manhattan(sys.argv[1])


'''
###############################################
#### Part four: Top gene list without SNPs ####
###############################################
'''
# start to log 
pre_message("Part four: Top gene list without SNPs ")

os.chdir(current_path)
sortTopGeneList(sys.argv[1])


'''
##########################################
## Part five: Top gene list with SNPs  ###
##########################################
'''
# start to log 
pre_message("Part five: Top gene list with SNPs ")

os.chdir(current_path)
sortTopGeneListWithSNPs(sys.argv[1])


'''
####################################
####  Part six:  Bubble plots   ####
####################################
'''
# start to log 
pre_message("Part six:  gene-tissue plots ")

os.chdir(current_path)
bubbleplot(sys.argv[1])


'''
####################################
#### Part seven: region plots   ####
####################################
'''

# start to log 
pre_message("Part seven: gene-level plots ")

os.chdir(current_path)
regionplot(sys.argv[1])


'''
######################################
#### Part eight: locuszoom plots  ####
######################################
'''
# start to log 
pre_message("Part eight: gwas snp region plots")

os.chdir(current_path)
locuszoom_plot(sys.argv[1])


'''
###################################
### write and close a log fle  ####
###################################
'''
# write logs 
write_logs()

# close logs 
finish_log()






