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
        if inputFilename == 'DGN-WB-unscaled.csv':
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

        if 'CrossTissue_elasticNet' in outputFileName:
            outputFileName = outputFileName[:-25]
        elif 'DGN-WB-unscaled' in outputFileName: 
            outputFileName = outputFileName [:-23]
        else: 
            outputFileName = outputFileName[:-25]
            outputFileName = outputFileName[3:]

        data.insert(2, 'tissue', outputFileName)
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

        if 'CrossTissue_elasticNet' in filename:
            tissue_name = filename[:-25]
        elif 'DGN-WB-unscaled' in filename: 
            tissue_name = filename [:-23]
        else: 
            tissue_name = filename[:-25]
            tissue_name = tissue_name[3:]
        df.insert(2, 'tissue', tissue_name)

        # sort data by defined column 
        df.sort_values(['pvalue', 'model_n'], ascending=[1, 0], inplace=True)
        total_rows = df.shape[0] # number of row count shape[1] is number of col count 
        total_rows_in_total += total_rows 
        cut_off = 0.05/total_rows 
        cut_off_list.append(cut_off)
        top_gene_list = df[df['pvalue'] < cut_off]
        # top_gene_list = top_gene_list[top_gene_list['pred_perf_R2'] > 0.01]
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
        # top_gene_list = top_gene_list[top_gene_list['pred_perf_R2'] > 0.01]

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


######################
###  Bubble Plot #####
######################

def bubbleplot(projectName):
    # load module 
    # load module 
    ggplot2 = importr('ggplot2')
    dplyr = importr('dplyr',  on_conflict="warn")

    currentPath = os.getcwd() 
    filePath_merged = currentPath + '/input/annotate/merged/'
    if not os.path.exists(filePath_merged): 
        msg = "Errors in finding path"
        print(msg)
    os.chdir(filePath_merged)

    data = r("data <- read.csv('mergedWithoutSNPs.csv')")
    robjects.globalenv['dataframe'] = data 

    # for snps data set 
    # set up file path 
    outputPath = currentPath + '/out/' + projectName + '/bubble-plots'
    if not os.path.exists(outputPath): 
         os.makedirs(outputPath)
    os.chdir(outputPath)

    # set up vectors for loop 
    startSites = r("startSites <- c(177043226, 156435952, 129541931, 82670771, 16914895, 136155000, 46500673, 43567337, 29180996, 17390291, 68021297)")
    snpsNames = r("snpsNames <- c('rs6755777', 'rs62274042', 'rs1400482', 'rs35094336', 'rs7032221', 'rs635634', 'rs7207826', 'rs1879586', 'rs62070645', 'rs4808075', 'DPEP2')")
    chrosome = r("chrosome <- c(2, 3, 8, 8, 9, 9, 17, 17, 17, 19, 16)")

    # for loop through each snps site +/- 1000000 bp  
    r("""
        # for loop through each snps site +/- 1000000 bp  
        for (i in 1:length(startSites)) 
        {
            message("TILE PLOTING ", i, ": ", snpsNames[i])

            # subset data for each snps 
            subData <- subset(data, data$start > startSites[i] - 1500000  & 
                data$start < startSites[i] + 1500000 & data$chr==chrosome[i])
            subData <- subData[order(subData$start),]
            subData <- mutate(subData, z_score=ifelse(subData$zscore > 0, '   +   ', '   -   ')) 
            # Labels = subData$gene_name
            subData$gene_name <- factor(subData$gene_name, levels=subData$gene_name)
            # subData$tissue <- factor(subData$tissue, levels=subData$tissue)


            # draw plot 
            p <- ggplot(subData, aes(x=subData$gene_name, y=subData$tissue, size=abs(subData$zscore/2)))
            p + 
            geom_point(aes(colour=z_score)) + 
            scale_color_manual(values=c('blue', 'brown')) +
            scale_size_continuous(guide=FALSE, range=c(0,max(abs(subData$zscore)/2)))+
            # ggtitle(paste('locus: ', snpsNames[i], '(chromosome', chrosome[i], ')')) +
            labs(x='Gene', y='Tissue') + 
            # scale_x_discrete(breaks = subData$gene_name, labels=Labels) + 
            theme(axis.text.x = element_text(size=10, face='bold', angle = 90, hjust = 1)) +
            theme(axis.text.y = element_text(size=10, face='bold')) +
            theme(plot.title = element_text(size=18, face='bold')) +
            # theme(axis.title= element_text(size=18, face='bold')) +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
            theme(legend.position = "none")

            # save plot 
            ggsave(paste('tile_', snpsNames[i], '.png', sep=''), width=10, height=10)
        } 
    """)

    print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

######################
###  Region Plot #####
######################

def regionplot(projectName):
    # load module 
    ggplot2 = importr('ggplot2')
    dplyr = importr('dplyr',  on_conflict="warn")

    # read merged data file 
        # read merged data file 
    currentPath = os.getcwd() 
    filePath_merged = currentPath + '/input/annotate/merged/'
    if not os.path.exists(filePath_merged): 
        msg = "Errors in finding path"
        add_log(msg)
    os.chdir(filePath_merged)

    data = r("data <- read.csv('mergedWithoutSNPs.csv')")
    robjects.globalenv['dataframe'] = data 

    # for snps data set 
    # set up file path 
    outputPath = currentPath + '/out/' + projectName + '/region-plots'
    if not os.path.exists(outputPath): 
         os.makedirs(outputPath)
    os.chdir(outputPath)

    # set up vectors for loop 
    startSites = r("startSites <- c(177043226, 156435952, 129541931, 82670771, 16914895, 136155000, 46500673, 43567337, 29180996, 17390291, 68021297)")
    snpsNames = r("snpsNames <- c('rs6755777', 'rs62274042', 'rs1400482', 'rs35094336', 'rs7032221', 'rs635634', 'rs7207826', 'rs1879586', 'rs62070645', 'rs4808075', 'DPEP2')")
    chrosome = r("chrosome <- c(2, 3, 8, 8, 9, 9, 17, 17, 17, 19, 16)")

    r("""
        # loop through snps site +/- 1000,000 bps 
        for (i in 1:length(startSites)) 
        {
        message("PLOTING ", i, ": ", snpsNames[i])

        # subset data 
        subData <- subset(data, data$start > startSites[i] - 1500000  & data$start < startSites[i] 
            + 1500000 & data$chr==chrosome[i])
        subData$logp <- -log10(subData$pvalue)
        subData <- subData[order(subData$start),]
        if (snpsNames[i] == 'rs635634') 
        {   # set up subset data 
            subData <- mutate(subData, sig=ifelse(subData$logp > -log10(0.05/nrow(data)), 'Most Sig', 
            ifelse(subData$logp > 5.30103 & subData$logp <= -log10(0.05/nrow(data)), 'Sig', 
            ifelse(subData$logp > -log10(0.05/nrow(subData)) 
                & subData$logp <= 5.30103, 'Less Sig','Not Sig')))) 
            subData$gene_name <- factor(subData$gene_name, levels=subData$gene_name)

            print(nrow(data))
            print(0.05/nrow(data))
            print(nrow(subData))
            print(0.05/nrow(subData))

            # draw plot
            p <- ggplot(subData, aes(x=gene_name, y=logp))
            p + geom_point(aes(colour = sig)) + 
            scale_color_manual(guide=FALSE, values=c('black', 'black', 'black', 'black')) +
                # ggtitle(paste("locus: ", snpsNames[i], '(chromosome', chrosome[i], ')')) + 
            labs(x='Gene', y='-log10(p-value)') +
            # geom_label_repel(Labels, aes(label=Labels)) +
            # geom_jitter(width = 0.5, height = 0.5) +
            geom_hline(yintercept = 5.30103, linetype='dashed', color='black') +   # 5 x 10-6 
            geom_hline(yintercept = -log10(0.05/nrow(data)), linetype='dashed', color='black') +
            geom_hline(yintercept = -log10(0.05/nrow(subData)), linetype='dashed', color='black') +
            # theme(axis.ticks = element_line(size = 0.1)) + 
            # theme(axis.ticks.length = unit(1, "cm")) + 
            theme(axis.text.x = element_text(size=12, face='bold', angle = 90, hjust = 1)) +
            theme(axis.text.y = element_text(size=12, face='bold')) + 
            theme(plot.title = element_text(size=18, face='bold')) +
            theme(axis.title.x = element_blank(), axis.title.y=element_text(size=18, face='bold')) + 
            theme(legend.position = "none")
        } 
        else if (snpsNames[i] == 'rs7032221'){
            subData <- mutate(subData, sig=ifelse(subData$logp > -log10(0.05/nrow(data)), 'Most Sig', 
            ifelse(subData$logp > 5.30103 & subData$logp <= -log10(0.05/nrow(data)), 'Sig', 
            ifelse(subData$logp > -log10(0.05/nrow(subData)) 
                & subData$logp <= 5.30103, 'Less Sig','Not Sig')))) 
            subData$gene_name <- factor(subData$gene_name, levels=subData$gene_name)

            print(nrow(data))
            print(0.05/nrow(data))
            print(nrow(subData))
            print(0.05/nrow(subData))

            p <- ggplot(subData, aes(x=gene_name, y=logp))
            p + geom_point(aes(colour = sig)) + 
            scale_color_manual(guide=FALSE, values=c('black', 'black', 'black', 'black')) +
                # ggtitle(paste("locus: ", snpsNames[i], '(chromosome', chrosome[i], ')')) + 
            labs(x='Gene', y='-log10(p-value)') +
            geom_hline(yintercept = 5.30103, linetype='dashed', color='black') +
            geom_hline(yintercept = -log10(0.05/nrow(data)), linetype='dashed', color='black') +
            geom_hline(yintercept = -log10(0.05/nrow(subData)), linetype='dashed', color='black') +
            theme(axis.text.x = element_text(size=12, face='bold', angle = 90, hjust = 1)) +
            theme(axis.text.y = element_text(size=12, face='bold')) + 
            theme(plot.title = element_text(size=18, face='bold')) +
            theme(axis.title.x = element_blank(), axis.title.y=element_text(size=18, face='bold')) + 
            theme(legend.position = "none")
        }
        else if (snpsNames[i] == 'rs1400482'){
            subData <- mutate(subData, sig=ifelse(subData$logp > -log10(0.05/nrow(data)), 'Most Sig', 
            ifelse(subData$logp > 5.30103 & subData$logp <= -log10(0.05/nrow(data)), 'Sig', 
            ifelse(subData$logp > -log10(0.05/nrow(subData)) 
                & subData$logp <= 5.30103, 'Less Sig','Not Sig')))) 
            subData$gene_name <- factor(subData$gene_name, levels=subData$gene_name)

            print(nrow(data))
            print(0.05/nrow(data))
            print(nrow(subData))
            print(0.05/nrow(subData))

            p <- ggplot(subData, aes(x=gene_name, y=logp))
            p + geom_point(aes(colour = sig)) + 
            scale_color_manual(guide=FALSE, values=c('black', 'black', 'black', 'black')) +
                # ggtitle(paste("locus: ", snpsNames[i], '(chromosome', chrosome[i], ')')) + 
            labs(x='Gene', y='-log10(p-value)') +
            geom_hline(yintercept = 5.30103, linetype='dashed', color='black') +
            geom_hline(yintercept = -log10(0.05/nrow(data)), linetype='dashed', color='black') +
            geom_hline(yintercept = -log10(0.05/nrow(subData)), linetype='dashed', color='black') +
            theme(axis.text.x = element_text(size=12, face='bold', angle = 90, hjust = 1)) +
            theme(axis.text.y = element_text(size=12, face='bold')) + 
            theme(plot.title = element_text(size=18, face='bold')) +
            theme(axis.title.x = element_blank(), axis.title.y=element_text(size=18, face='bold')) + 
            theme(legend.position = "none")
        }
        else {
            subData <- mutate(subData, sig=ifelse(subData$logp > -log10(0.05/nrow(data)), 'Most Sig', 
            ifelse(subData$logp > 5.30103 & subData$logp <= -log10(0.05/nrow(data)), 'Sig', 
            ifelse(subData$logp > -log10(0.05/nrow(subData)) 
                & subData$logp <= 5.30103, 'Less Sig','Not Sig')))) 
            subData$gene_name <- factor(subData$gene_name, levels=subData$gene_name)

            print(nrow(data))
            print(0.05/nrow(data))
            print(nrow(subData))
            print(0.05/nrow(subData))

            p <- ggplot(subData, aes(x=gene_name, y=logp))
            p + geom_point(aes(colour = sig)) + 
            scale_color_manual(guide=FALSE, values=c('black', 'black', 'black', 'black')) +
                # ggtitle(paste("locus: ", snpsNames[i], '(chromosome', chrosome[i], ')')) + 
            labs(x='Gene', y='-log10(p-value)') +
            geom_hline(yintercept = 5.30103, linetype='dashed', color='black') +
            geom_hline(yintercept = -log10(0.05/nrow(data)), linetype='dashed', color='black') +
            geom_hline(yintercept = -log10(0.05/nrow(subData)), linetype='dashed', color='black') +
            theme(axis.text.x = element_text(size=12, face='bold', angle = 90, hjust = 1)) +
            theme(axis.text.y = element_text(size=12, face='bold')) + 
            theme(plot.title = element_text(size=18, face='bold')) +
            theme(axis.title.x = element_blank(), axis.title.y=element_text(size=18, face='bold')) + 
            theme(legend.position = "none")
        } 

        # output plot 
        ggsave(paste(snpsNames[i], '.png',sep=''), width=10, height=10)
        } 
    """) 

    print(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")


