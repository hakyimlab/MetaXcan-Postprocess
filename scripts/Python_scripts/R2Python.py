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

total_cut_off = 0 
cut_off_list = [] 

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
        if inputFilename == 'DGN-WB-unscale.csv':
            annotatedData = r("inner_join(data[, c('gene', 'gene_name', 'zscore', 'pvalue', 'model_n', 'pred_perf_R2')], grch37[, c('symbol', 'chr', 'start', 'end')], by=c('gene'='symbol'))")
        else:
            annotatedData = r("inner_join(data[, c('gene', 'gene_name', 'zscore', 'pvalue', 'model_n', 'pred_perf_R2')], grch37[, c('ensgene', 'chr', 'start', 'end')], by =c('gene'='ensgene'))")
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

    # for snps data set 
    # set up file path 
    currentPath = get_current_path()
    outputPath = os.getcwd() + '/out/' + projectName + '/QQ-Plot'
    if not os.path.exists(outputPath): 
         os.makedirs(outputPath)
    filePath = os.getcwd() + '/out/' + projectName 
    os.chdir(filePath)


    # get file lists 
    outputFileList = glob.glob("*.csv")

    dfList8 = []
    for outputFileName in outputFileList: 
        # get GWAS output data 
        print("CONCATING FILE: " + outputFileName)
        r("data <- read.csv('%s')" %outputFileName)
        data = r('data <- na.omit(data)')
        robjects.globalenv['dataframe'] = data
        dfList8.append(data)
    concatDf8 = pandas.concat(dfList8, axis = 0)

    #output data 
    mergedOutputPath = currentPath + '/out/' + projectName + '/merged/'
    if not os.path.exists(mergedOutputPath): 
         os.makedirs(mergedOutputPath)
    concatDf8.to_csv(mergedOutputPath + 'merged.csv', index = None)
    print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")


    
    # loop through file lists 
    for outputFileName in outputFileList: 
        print("QQ-PLOT(snps): " + outputFileName)

        # get GWAS output data 
        r("data <- read.csv('%s')" %outputFileName)
        data = r('data <- na.omit(data)')
        robjects.globalenv['dataframe'] = data

        os.chdir(outputPath)

        # draw qq-plot and save them to files 
        r.png('%s%s%s'%('QQ-Plot_', outputFileName[:-4],'.png'))
        qqman.qq(data['pvalue'],  main = "Q-Q plot of GWAS p-values")
        r['dev.off']()
        os.chdir(filePath)

        print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

    print("QQ-PLOT(snps): " + "merged.csv")

    # get GWAS output data 
    os.chdir(mergedOutputPath)
    r("data8 <- read.csv('merged.csv')")
    data8 = r('data8 <- na.omit(data)')
    robjects.globalenv['dataframe'] = data8

    os.chdir(outputPath)

    # draw qq-plot and save them to files 
    r.png('%s%s%s'%('QQ-Plot_', outputFileName[:-4],'.png'))
    qqman.qq(data['pvalue'],  main = "Q-Q plot of GWAS p-values")
    r['dev.off']()
    os.chdir(filePath)

    print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

    # for without snps data set 
    # set up file path 
    outputPathWithoutSNPs = currentPath + '/input/annotate' +'/QQ-Plot'
    if not os.path.exists(outputPathWithoutSNPs): 
         os.makedirs(outputPathWithoutSNPs)
    filePathWithoutSNPs = currentPath + '/input/annotate/'
    os.chdir(filePathWithoutSNPs)

    # get file lists 
    outputFileListWithoutSNPs = glob.glob("*.csv")

    dfListWithoutSNPs = []
    for outputFileName in outputFileListWithoutSNPs: 
        # get GWAS output data 
        print("CONCATING FILE: " + outputFileName)
        r("data <- read.csv('%s')" %outputFileName)
        data = r('data <- na.omit(data)')
        robjects.globalenv['dataframe'] = data
        data.insert(2, 'tissue', outputFileName[:-25])
        dfListWithoutSNPs.append(data)
    concatDfWithoutSNPs = pandas.concat(dfListWithoutSNPs, axis = 0)

    #output data 
    mergedpath = filePathWithoutSNPs + 'merged/'
    if not os.path.exists(mergedpath): 
        os.makedirs(mergedpath)
    concatDfWithoutSNPs.to_csv(mergedpath + 'mergedWithoutSNPs.csv', index = None)
    print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

    # loop through file lists 
    for outputFileName in outputFileListWithoutSNPs: 
        print("QQ-PLOT(no snps): " + outputFileName)

        # get GWAS output data 
        r("data <- read.csv('%s')" %outputFileName)
        data = r('data <- na.omit(data)')
        robjects.globalenv['dataframe'] = data
        os.chdir(outputPathWithoutSNPs)

        # draw qq-plot and save them to files 
        r.png('%s%s%s'%('QQ-Plot_', outputFileName[:-4],'.png'))
        qqman.qq(data['pvalue'],  main = "Q-Q plot of GWAS p-values")
        r['dev.off']()
        os.chdir(filePathWithoutSNPs)

        print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

    print("QQ-PLOT(no snps): " + "mergedWithoutSNPs.csv")

    # get GWAS output data 
    os.chdir(filePathWithoutSNPs + 'merged/')
    r("data <- read.csv('mergedWithoutSNPs.csv')")
    data = r('data <- na.omit(data)')
    robjects.globalenv['dataframe'] = data
    os.chdir(outputPathWithoutSNPs)

    # draw qq-plot and save them to files 
    r.png('%s%s%s'%('QQ-Plot_', 'mergedWithoutSNPs', '.png'))
    qqman.qq(data['pvalue'],  main = "Q-Q plot of GWAS p-values")
    r['dev.off']()
    os.chdir(filePathWithoutSNPs)

    print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")


###########################
####  Manhattan-Plot  #####
########################### 

def manhattan(projectName):
    # load module 
    qqman = importr('qqman')

    currentPath = get_current_path()

    # for snps data set 
    # set up file path 
    outputPath = os.getcwd() + '/out/' + projectName +'/Manhattan-Plot'

    if not os.path.exists(outputPath): 
         os.makedirs(outputPath)

    filePath = os.getcwd() + '/out/' + projectName 
    os.chdir(filePath)

    # get file lists 
    outputFileList = glob.glob("*.csv")

    # loop through file lists 
    for outputFileName in outputFileList: 
        print("MANHATTAN-PLOT(snps): " + outputFileName)

        # get GWAS output data 
        r("data <- read.csv('%s')" %outputFileName)
        r('data$chr <- as.numeric(as.character(data$chr))')
        data = r('data <- na.omit(data)')
        robjects.globalenv['dataframe'] = data

        os.chdir(outputPath)

        # draw manhattan and save them to files 
        r.png('%s%s%s'%('Manhattan-Plot_', outputFileName[:-4], '.png'))
        qqman.manhattan(data, chr = 'chr', bp='start', p='pvalue', snp='rsid', main = 'Manhattan plot')
        r['dev.off']()
        os.chdir(filePath)

        print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

    print("MANHATTAN-PLOT(snps): " + 'merged.csv')

    # get GWAS output data 
    os.chdir(filePath + '/merged/')
    r("data8 <- read.csv('merged.csv')")
    r('data8$chr <- as.numeric(as.character(data8$chr))')
    data8 = r('data8 <- na.omit(data8)')
    robjects.globalenv['dataframe'] = data8

    os.chdir(outputPath)

    # draw manhattan and save them to files 
    r.png('%s%s%s'%('Manhattan-Plot_', 'merged.csv', '.png'))
    qqman.manhattan(data, chr = 'chr', bp='start', p='pvalue', snp='rsid', main = 'Manhattan plot')
    r['dev.off']()
    os.chdir(filePath)

    print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")


    # for raw data set without snps  
    outputPathWithoutSNPs = currentPath + '/input/annotate' +'/Manhattan-Plot'

    if not os.path.exists(outputPathWithoutSNPs): 
         os.makedirs(outputPathWithoutSNPs)
    filePathWithoutSNPs = currentPath + '/input/annotate/'
    os.chdir(filePathWithoutSNPs)

    # get file lists 
    outputFileListWithoutSNPs = glob.glob("*.csv")

    # loop through file lists 
    for outputFileName in outputFileListWithoutSNPs: 
        print("MANHATTAN-PLOT(no snps): " + outputFileName)

        # get GWAS output data 
        r("data <- read.csv('%s')" %outputFileName)
        r('data$chr <- as.numeric(as.character(data$chr))')
        data = r('data <- na.omit(data)')
        robjects.globalenv['dataframe'] = data

        os.chdir(outputPathWithoutSNPs)

        # draw manhattan and save them to files 
        r.png('%s%s%s'%('Manhattan-Plot_', outputFileName[:-4], '.png'))
        qqman.manhattan(data, chr = 'chr', bp='start', p='pvalue', snp='rsid', main = 'Manhattan plot')
        r['dev.off']()
        os.chdir(filePathWithoutSNPs)

        print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

    print("MANHATTAN-PLOT(no snps): " + "mergedWithoutSNPs.csv")

    # get GWAS output data 
    os.chdir(filePathWithoutSNPs + 'merged/')
    r("data <- read.csv('mergedWithoutSNPs.csv')")
    r('data$chr <- as.numeric(as.character(data$chr))')
    data = r('data <- na.omit(data)')
    robjects.globalenv['dataframe'] = data

    os.chdir(outputPathWithoutSNPs)

    # draw manhattan and save them to files 
    r.png('%s%s%s'%('Manhattan-Plot_', 'mergedWithoutSNPs.csv', '.png'))
    qqman.manhattan(data, chr = 'chr', bp='start', p='pvalue', snp='rsid', main = 'Manhattan plot')
    r['dev.off']()
    os.chdir(filePathWithoutSNPs)

    print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

######################################
####   Top Gene List without SNPs ####
###################################### 

def sortTopGeneList(projectName):

    currentPath = get_current_path()

    os.chdir(currentPath + '/input/annotate/')

    fileList = glob.glob("*.csv")

    out_path = currentPath + '/input/annotate/sorted/'
    if not os.path.exists(out_path): 
        os.makedirs(out_path)

    dfList = []
    total_rows_in_total = 0 
    for filename in fileList:
        print("TOP GENE LIST (no snps): " + filename)
        df = pandas.read_csv(filename) 
        df.insert(2, 'tissue', filename[:-25])

        # sort data by defined column 
        df.sort_values(['pvalue', 'model_n'], ascending=[1, 0], inplace=True)
        total_rows = df.shape[0] # number of row count shape[1] is number of col count 
        total_rows_in_total += total_rows 
        cut_off = 0.05/total_rows 
        cut_off_list.append(cut_off)
        top_gene_list = df[df['pvalue'] < cut_off]
        top_gene_list = top_gene_list[top_gene_list['pred_perf_R2'] > 0.01]
        # top_gene_list.drop('gene', axis=1, inplace=True)
        print('CALCULATING cut_off p-Value: ' + str(cut_off))

        #output data 
        top_gene_list.to_csv(out_path + filename, index = None)

        dfList.append(top_gene_list)

        print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

    print("TOP GENE LIST(all without snps): " + "sortedTopGenesAllTissues.csv")
    concatDf = pandas.concat(dfList, axis = 0)

    global total_cut_off 
    total_cut_off = 0.05/total_rows_in_total 
    print('CALCULATING total cut_off p-Value: ' + str(total_cut_off))

    concatDf = concatDf[concatDf['pvalue'] < total_cut_off]

    # sort data by defined column 
    concatDf.sort_values(['pvalue', 'model_n'], ascending=[1, 0], inplace=True)

    #output data 
    concatDf.to_csv(out_path + 'sortedTopGenesAllTissues.csv', index = None)


    print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")


####################################
####   Top Gene List with SNPs #####
#################################### 

def sortTopGeneListWithSNPs(projectName):

    currentPath = get_current_path()

    os.chdir(currentPath + '/out/' + projectName)

    fileList = glob.glob("*.csv")

    out_path = currentPath + '/out/' + projectName + '/sorted/'
    if not os.path.exists(out_path): 
        os.makedirs(out_path)

    dfList = []
    for i in range(len(fileList)):
        print("TOP GENE LIST (snps): " + fileList[i])
        df = pandas.read_csv(fileList[i]) 

        # sort data by defined column 
        df.sort_values(['pvalue', 'model_n'], ascending=[1, 0], inplace=True)
        # total_rows = df.shape[0] # number of row count shape[1] is number of col count 
        cut_off = cut_off_list[i]
        top_gene_list = df[df['pvalue'] < cut_off]
        top_gene_list = top_gene_list[top_gene_list['pred_perf_R2'] > 0.01]

        # top_gene_list.drop('gene', axis=1, inplace=True)
        print('CALCULATING cut_off p-Value: ' + str(cut_off))

        #output data 
        top_gene_list.to_csv(out_path + fileList[i], index = None)

        dfList.append(top_gene_list)

        print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

    print("TOP GENE LIST(all with snps): " + "sortedTopGenesAllTissues.csv")
    concatDf = pandas.concat(dfList, axis = 0)
    global total_cut_off
    print('Total_cut_off p-Value: ' + str(total_cut_off))
    concatDf = concatDf[concatDf['pvalue'] < total_cut_off]

    # sort data by defined column 
    concatDf.sort_values(['pvalue', 'model_n'], ascending=[1, 0], inplace=True)

    #output data 
    concatDf.to_csv(out_path + 'sortedTopGenesAllTissues.csv', index = None)


    print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")








