''' This script wraps R packages such as annotables and qqman in python environment using rpy2 '''
''' Please see more documentation at http://rpy2.readthedocs.org/en/version_2.7.x/ ''' 


################################
#### Software prerequistes #####
################################

# rpy2: pip install rpy2
# annotables: install.packages('annotables')
# dplyr: install.packages('dplyr')

###################
#### Modules  #####
###################

from rpy2.robjects import r
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2 import robjects
from helpers import * 
from datetime import datetime 
import os, pandas, glob

pandas2ri.activate()


######################
#### Annotations #####
######################

def annotate_metaxcan_result():
    # load modules 
    importr('annotables')
    importr('dplyr',  on_conflict="warn")

    # setup path 
    currentPath = get_current_path()
    inputpath = get_input_path_without_filename()

    annoateInputPath = inputpath + '/annotate/'
    if not os.path.exists(annoateInputPath): 
        os.makedirs(annoateInputPath)

    os.chdir(inputpath)

    # get file lists 
    inputFileList = glob.glob("*.csv")

    # loop through file lists 
    for inputFilename in inputFileList:
        print("ANNOTATING GENES: " + inputFilename)

        # GWAS data 
        r("data <- read.csv('%s')" %inputFilename)
        data = r('na.omit(data)')
        robjects.globalenv['dataframe'] = data

        # annotation data 
        grch37 = r('grch37')
        robjects.globalenv['dataframe'] = grch37

        # annotate GWAS data 
        annotatedData = r("inner_join(data[, c('gene', 'gene_name', 'zscore', 'pvalue', 'model_n')], grch37[, c('ensgene', 'chr', 'start', 'end')], by =c('gene'='ensgene'))")
        annotatedData = annotatedData.drop_duplicates()

        # ouput annotated data 
        annotatedData.to_csv(annoateInputPath + inputFilename[:-4] + '_annotated.csv', index=None)

        print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

    os.chdir(annoateInputPath)
    

####################
####  QQ-Plot  #####
####################

def qqplot(projectName):
    # load module 
    qqman = importr('qqman')

    # set up file path 
    currentPath = get_current_path()
    outputPath = os.getcwd() + '/out/' + projectName + '/QQ-Plot'
    if not os.path.exists(outputPath): 
         os.makedirs(outputPath)
    filePath = os.getcwd() + '/out/' + projectName 
    os.chdir(filePath)

    # get file lists 
    outputFileList = glob.glob("*.csv")

    # loop through file lists 
    for outputFileName in outputFileList: 
        print("QQ-PLOT: " + outputFileName)

        # get GWAS output data 
        r("data <- read.csv('%s')" %outputFileName)
        data = r('na.omit(data)')
        robjects.globalenv['dataframe'] = data
        os.chdir(outputPath)

        # draw qq-plot and save them to files 
        r.png('%s%s'%(outputFileName[:-4],'.png'), width=300, height=300)
        qqman.qq(data['pvalue'],  main = "Q-Q plot of GWAS p-values")
        r['dev.off']()
        os.chdir(filePath)

        print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")


###########################
####  Manhattan-Plot  #####
########################### 

def manhattan(projectName):
    # load module 
    qqman = importr('qqman')

    # set up file path 
    currentPath = get_current_path()
    outputPath = os.getcwd() + '/out/' + projectName +'/Manhattan-Plot'

    if not os.path.exists(outputPath): 
         os.makedirs(outputPath)

    filePath = os.getcwd() + '/out/' + projectName 
    os.chdir(filePath)

    # get file lists 
    outputFileList = glob.glob("*.csv")

    # loop through file lists 
    for outputFileName in outputFileList: 
        print("MANHATTAN-PLOT: " + outputFileName)

        # get GWAS output data 
        r("data <- read.csv('%s')" %outputFileName)
        data = r('na.omit(data)')
        robjects.globalenv['dataframe'] = data

        os.chdir(outputPath)

        # draw manhattan and save them to files 
        r.png('%s%s'%(outputFileName[:-4], '.png'), width=300, height=300)
        qqman.manhattan(data, chr = 'chr', bp='start', p='pvalue', snp='rsid', main = 'Manhattan plot')
        r['dev.off']()
        os.chdir(filePath)

        print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")




