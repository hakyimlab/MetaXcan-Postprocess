# A simple "main()" that loads MetaXcanPostprocessing file 

import sys, os, glob 
from helpers import *
from MetaXcanPostprocessing import * 
from R2Python import * 

if len(sys.argv) != 2:
    hints1 = "Syntax: "
    hints2 = "    python run_MetaXcanPostprocessing.py <project name>"
    exit()

current_path = os.getcwd()

# Run metaxcan postprocessing 

#############################
### Part one: Annotation ####
#############################

# Annotated files
print (LINE)
print ('      Getting started MetaXcan postprocessing -- Part one: Annotating genes ')
print (LINE)
annotate_metaxcan_result()


#################################
#### Part two: Fetching SNPs ####
#################################

print (LINE)
print ('      Getting started MetaXcan postprocessing -- Part two: Fetching gene SNPs  ')
print (LINE)

inputFileList = glob.glob("*.csv")

os.chdir(current_path)

for inputFileName in inputFileList: 
    inputName = inputFileName[:-4]
    fullInputName = sys.argv[1] + '_' + inputName 
    
    # Start up the logfile
    open_log(fullInputName)

    # Init class and read .csv file 
    snps = TopGeneListSnps(inputName, fullInputName) 

    # Fetch top gene list from sqlite databases 
    snps.fetchTopGeneList(inputName) 

    # Output data
    snps.outputTopGeneList(fullInputName) 

     # Log data 
    snps.logTopGeneList(inputName)

    # Close the logfile
    finish_log()

    print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")


#############################
#### Part three: QQ-Plot ####
#############################

print (LINE)
print ('      Getting started MetaXcan postprocessing -- Part three: QQ-Plot  ')
print (LINE)

os.chdir(current_path)
qqplot(sys.argv[1]) 


###################################
#### Part four: Manhattan-Plot ####
###################################

print (LINE)
print ('      Getting started MetaXcan postprocessing -- Part four: Manhattan-Plot  ')
print (LINE)

os.chdir(current_path)
manhattan(sys.argv[1])


###############################################
#### Part five: Top gene list without SNPs ####
###############################################

print (LINE)
print ('      Getting started MetaXcan postprocessing -- Part five: Top gene list without SNPs ')
print (LINE)

os.chdir(current_path)
sortTopGeneList(sys.argv[1])


############################################
#### Part Six: Top gene list with SNPs  ####
############################################

print (LINE)
print ('      Getting started MetaXcan postprocessing -- Part six: Top gene list with SNPs ')
print (LINE)

os.chdir(current_path)
sortTopGeneListWithSNPs(sys.argv[1])

####################################
#### Part Seven: Bubble plots   ####
####################################

print (LINE)
print ('      Getting started MetaXcan postprocessing -- Part seven: Bubble plot ')
print (LINE)

os.chdir(current_path)
bubbleplot(sys.argv[1])

####################################
#### Part Eight: region plots   ####
####################################

print (LINE)
print ('      Getting started MetaXcan postprocessing -- Part eight: Region plot ')
print (LINE)

os.chdir(current_path)
regionplot(sys.argv[1])


