''' 
This python script wraps some R packages such as annotables and qqman, or some R codes in python environment using rpy2
Please see more documentation at http://rpy2.readthedocs.org/en/version_2.7.x/ 
''' 

'''
################################
#### Software prerequistes #####
################################
'''
# python packages: 
# rpy2:       pip install rpy2

# R packages: 
# annotables:     install.packages("devtools")
#                 devtools::install_github("stephenturner/annotables")
# dplyr:          install.packages('dplyr')
# ggplot2:        install.packages("tidyverse") 
# qqman:          install.packages("qqman") 


'''
###################
#### Modules  #####
###################
'''
from rpy2.robjects import r
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2 import robjects
from helpers import * 
from datetime import datetime 
import os, pandas, glob, sqlite3, shutil, subprocess
import uuid as myuuid 

pandas2ri.activate()

projectID = str(myuuid.uuid4())


'''
######################
#### Annotations #####
######################
'''
def annotate_metaxcan_result(projectName):

    global projectID 

    # load R packages  
    importr('annotables')
    importr('dplyr',  on_conflict="warn")

    # current path 
    currentPath = os.getcwd()

    # root path 
    os.chdir('../../../')
    root_path = os.getcwd()

    # input path 
    input_path = root_path + '/input/'

    if not os.path.exists(input_path): 
        warning = "Please make sure that you have created an input folder with input files"
        add_log(warning)
        add_log('Input path should be: %s' %input_path)
    os.chdir(input_path)

    # annoated_files output file path
    annoated_files_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME + '/annotated_metaxcan_output_files/'
    if not os.path.exists(annoated_files_path): 
        os.makedirs(annoated_files_path)

    # get input file lists (raw metaxcan output files *.csv) 
    inputFileList = glob.glob("*.csv")

    # loop through file lists 
    for inputFilename in inputFileList:
        msg = "ANNOTATING GENES: " + inputFilename
        add_log(msg)

        # read data 
        r("data <- read.csv('%s')" %inputFilename)
        data = r('na.omit(data)')
        robjects.globalenv['dataframe'] = data

        # annotation library 
        grch37 = r('grch37')
        robjects.globalenv['dataframe'] = grch37

        # annotating  
        # if inputFilename == 'xxxxxDGN-WB-unscaled.csv':
        #     annotatedData = r("inner_join(data[, c('gene', 'gene_name', 'zscore', 'effect_size', 'pvalue', 'n_snps_in_model', 'pred_perf_p', 'pred_perf_R2')], grch37[, c('symbol', 'chr', 'start', 'end')], by=c('gene'='symbol'))")
        # else:
        annotatedData = r("inner_join(data[, c('gene', 'gene_name', 'zscore', 'effect_size', 'pvalue', 'n_snps_in_model', 'pred_perf_p', 'pred_perf_R2')], grch37[, c('ensgene', 'chr', 'start', 'end')], by =c('gene'='ensgene'))")
        annotatedData = annotatedData.drop_duplicates()

        # ouput annotated data 
        os.chdir(annoated_files_path)
        annotatedData.to_csv(annoated_files_path + inputFilename[:-4] + '_annotated.csv', index=None)

        os.chdir(input_path)

    # merge all annotated metaxcan files 
    # get annotated file lists 
    os.chdir(annoated_files_path)
    outputFileListWithoutSNPs = glob.glob("*.csv")

    # merged annotated file path 
    merged_annotated_file_path = annoated_files_path + 'merged/'
    if not os.path.exists(merged_annotated_file_path): 
        os.makedirs(merged_annotated_file_path)

    dfListWithoutSNPs = []
    for outputFileName in outputFileListWithoutSNPs: 
        # read data 
        msg = "MERGING ALL ANNOTATED FILES ..."
        add_log(msg)
        msg = "CONCATING FILE: " + outputFileName
        add_log(msg)
        
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
    os.chdir(merged_annotated_file_path)
    concatDfWithoutSNPs.to_csv(merged_annotated_file_path + 'merged_annotated_metaxcan_output.csv', index = None)

    add_log(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")
    

'''
####################
####  QQ-Plot  #####
####################
'''
def qqplot(projectName):

    global projectID 

    # load R package
    qqman = importr('qqman')
    
    # current path 
    currentPath = os.getcwd()

    # back to root path 
    os.chdir('../../../')
    root_path = os.getcwd()

    # annotated metaxcan file path 
    annoated_files_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME + '/annotated_metaxcan_output_files/'
    os.chdir(annoated_files_path)

    qqplot_output_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME +'/QQ-plot/'
    if not os.path.exists(qqplot_output_path): 
         os.makedirs(qqplot_output_path)

    outputFileListWithoutSNPs = glob.glob("*.csv")

    # loop through file lists 
    for outputFileName in outputFileListWithoutSNPs: 
        msg = "QQ-PLOT(no snps): " + outputFileName
        add_log(msg)

        # get GWAS output data 
        r("data <- read.csv('%s')" %outputFileName)
        data = r('data <- na.omit(data)')
        robjects.globalenv['dataframe'] = data
        os.chdir(qqplot_output_path)

        # draw qq-plot and save them to files 
        r.png('%s%s%s'%('QQ-Plot_', outputFileName[:-4],'.png'))
        qqman.qq(data['pvalue'],  main = "MetaXcan genes")
        r['dev.off']()
        os.chdir(annoated_files_path)

    msg = "QQ-PLOT(no snps): " + "merged_annotated_metaxcan_output"
    add_log(msg)

    # set up path and read merged annotated metaxcan output files  
    os.chdir(annoated_files_path + 'merged/')
    r("data <- read.csv('merged_annotated_metaxcan_output.csv')")
    data = r('data <- na.omit(data)')
    robjects.globalenv['dataframe'] = data
    os.chdir(qqplot_output_path)

    # draw qq-plot and save them to files 
    r.png('%s%s%s'%('QQ-Plot_', 'metaxcan_output', '.png'))
    qqman.qq(data['pvalue'],  main = "QQ plot")
    r['dev.off']()

    add_log(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")


'''
###########################
####  Manhattan-Plot  #####
########################### 
'''
def manhattan(projectName):

    global projectID 

    # load module 
    qqman = importr('qqman')

    # current path 
    currentPath = os.getcwd()

    # back to root path 
    os.chdir('../../../')
    root_path = os.getcwd()

    # annotated metaxcan file path 
    annoated_files_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME +'/annotated_metaxcan_output_files/'
    os.chdir(annoated_files_path)

    manhattanplot_output_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME + '/Manhattan-Plot/'
    if not os.path.exists(manhattanplot_output_path): 
         os.makedirs(manhattanplot_output_path)

    # get file lists 
    outputFileListWithoutSNPs = glob.glob("*.csv")

    # loop through file lists 
    for outputFileName in outputFileListWithoutSNPs: 
        msg = "MANHATTAN-PLOT(no snps): " + outputFileName
        add_log(msg)

        # read data 
        r("data <- read.csv('%s')" %outputFileName)
        r('data$chr <- as.numeric(as.character(data$chr))')
        data = r('data <- na.omit(data)')
        robjects.globalenv['dataframe'] = data

        os.chdir(manhattanplot_output_path)

        # draw manhattan and save them to files 
        r.png('%s%s%s'%('Manhattan-Plot_', outputFileName[:-4], '.png'))
        qqman.manhattan(data, chr = 'chr', bp='start', p='pvalue', snp='gene_name', main = 'MetaXcan genes')
        r['dev.off']()
        os.chdir(annoated_files_path)

    msg = "MANHATTAN-PLOT(no snps): " + "merged_annotated_metaxcan_output"
    add_log(msg)

    # read data 
    os.chdir(annoated_files_path + 'merged/')
    r("data <- read.csv('merged_annotated_metaxcan_output.csv')")
    r('data$chr <- as.numeric(as.character(data$chr))')
    data = r('data <- na.omit(data)')
    robjects.globalenv['dataframe'] = data

    os.chdir(manhattanplot_output_path)

    # draw manhattan and save them to files 
    r.png('%s%s%s'%('Manhattan-Plot_', 'merged_annotated_metaxcan_output.csv', '.png'))
    qqman.manhattan(data, chr = 'chr', bp='start', p='pvalue', snp='gene_name', main = 'Manhattan plot')
    r['dev.off']()

    add_log(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

'''
######################################
####   Top Gene List without SNPs ####
###################################### 
'''
def sortTopGeneList(projectName):

    global projectID 

    # current path 
    currentPath = os.getcwd()

    # back to root path 
    os.chdir('../../../')
    root_path = os.getcwd()

    # annotated metaxcan file path 
    annoated_files_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME +'/annotated_metaxcan_output_files/'
    os.chdir(annoated_files_path)

    top_genes_output_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME +'/top_genes/'
    if not os.path.exists(top_genes_output_path): 
         os.makedirs(top_genes_output_path)

    # get file lists 
    fileLists = glob.glob("*.csv")

    # dfList = []
    # total_rows_in_total = 0 
    for filename in fileLists:
        msg = "TOP GENE LIST: " + filename
        add_log(msg)

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
        df.sort_values(['pvalue', 'n_snps_in_model'], ascending=[1, 0], inplace=True)
        total_rows = df.shape[0] # number of row count shape[1] is number of col count 
        msg = 'total_row within tissue: ' + str(total_rows)
        add_log(msg)
        # total_rows_in_total += total_rows 
        cut_off = 0.05/total_rows 

        # cut_off_list.append(cut_off)
        top_gene_list = df[df['pvalue'] < cut_off]
        # top_gene_list = top_gene_list[top_gene_list['pred_perf_R2'] > 0.01]
        # top_gene_list.drop('gene', axis=1, inplace=True)
        msg = 'CALCULATING TISSUE-WIDE P-VALUE: ' + str(cut_off)
        add_log(msg)
        add_log(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")

        #output data 
        top_gene_list.to_csv(top_genes_output_path + "sorted_top_genes_%s"%filename, index = None)

        # dfList.append(top_gene_list)

    os.chdir(annoated_files_path + 'merged/')

    msg = "TOP GENES: merged_annotated_metaxcan_output.csv"
    add_log(msg)

    df = pandas.read_csv('merged_annotated_metaxcan_output.csv') 

    total_rows_in_total = df.shape[0] 
    msg = 'total_row_in_total: ' + str(total_rows_in_total)
    add_log(msg)

    total_cut_off = 0.05/total_rows_in_total 

    msg = 'CALCULATING GENOME-WIDE P-VALUE: ' + str(total_cut_off)
    add_log(msg)

    df = df[df['pvalue'] < total_cut_off]

    # sort data by defined column 
    df.sort_values(['pvalue', 'n_snps_in_model'], ascending=[1, 0], inplace=True)

    #output data 
    df.to_csv(top_genes_output_path+ 'sorted_top_genes.csv', index = None)

    add_log(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")


'''
####################################
####   Top Gene List with SNPs #####
#################################### 
'''
def sortTopGeneListWithSNPs(projectName):

    global projectID 

    # current path 
    currentPath = os.getcwd()

    # back to root path 
    os.chdir('../../../')
    root_path = os.getcwd()

    # set up path 
    top_genes_file_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME + '/top_genes/'
    os.chdir(top_genes_file_path)

    # read data 
    top_genes = pandas.read_csv('sorted_top_genes.csv')

    # add into files 
    gene_lists = top_genes['gene_name']
    tissue_lists = top_genes['tissue']
    pvalue_lists = top_genes['pvalue']
    zscore_lists = top_genes['zscore']
    model_n_lists = top_genes['n_snps_in_model']
    pred_perf_R2_lists = top_genes['pred_perf_R2']
    chr_lists = top_genes['chr']
    start_lists = top_genes['start']
    end_lists = top_genes['end']

    # get db lists 
    input_path = root_path + '/input/'

    if not os.path.exists(input_path): 
        warning = "Please make sure that you have created an input folder with input files"
        add_log(warning)
        add_log('Input path should be: %s' %input_path)
    os.chdir(input_path)

    dbFileList = glob.glob("*.db")

    database_names = []
    for dbFilename in dbFileList:
       database_names.append(dbFilename)

    # Loop through databases 
    query_output_list = []
    for i in range(len(database_names)):
        for k in range(len(tissue_lists)): 
            if tissue_lists[k] in database_names[i]: 
                # Connect databases 
                conn = sqlite3.connect(database_names[i]) 

                full_query_name = None 
                if database_names[i] == 'DGN-WB-unscaled_0.5.db':
                    full_query_name = SQL_QUERY_PREFIX_DNG + gene_lists[k] + "'"
                else: 
                    full_query_name = SQL_QUERY_PREFIX + gene_lists[k] + "'"

                query_output = pandas.read_sql(full_query_name, conn, index_col=None)

                if database_names[i] == 'DGN-WB-unscaled_0.5.db':
                    query_output.rename(columns={'gene':'genename'}, inplace=True) 

                # # Add correspinding parameters to the new output file 
                query_output['tissue'] = tissue_lists[k]
                query_output['pvalue'] =  pvalue_lists[k]
                query_output['zscore'] =  zscore_lists[k]
                query_output['model_n'] = model_n_lists[k]
                query_output['pred_perf_R2'] = pred_perf_R2_lists[k]
                query_output['chr'] =  chr_lists[k]
                query_output['start'] = start_lists[k]
                query_output['end'] = end_lists[k]
                query_output_list.append(query_output) 

                msg = 'FETCH SNPs for GENE %s ' % gene_lists[k] +  'FROM DATABASE: %s' % database_names[i]
                add_log(msg)

                # Close database
                conn.close() 

    # Merge output data 
    query_output_of = pandas.concat(query_output_list, axis = 0)   
    # gwas_snp = pandas.read_csv('gwas_snp.txt', sep="\s+")
    # gwas_snp.rename(columns={'MarkerName':'rsid', 'P-value':'gwas_pvalue'}, inplace=True)

    # query_output_of.merge(gwas_snp, on='rsid', how='inner')

    # get tissue abbr name 
    # load module 
    ggplot2 = importr('ggplot2')
    dplyr = importr('dplyr',  on_conflict="warn")

    r("""

    tissue_abbr <- read.delim('gtex_tissue_abbr.txt', sep='\t') %>% 
    select(tissue_site_detail_abbr,tissue_site_detail_id)
    colnames(tissue_abbr)[colnames(tissue_abbr) == 'tissue_site_detail_id'] = "tissue"

    """) 


    top_genes_snps_output_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME + '/top_genes_snps/'

    if not os.path.exists(top_genes_snps_output_path): 
         os.makedirs(top_genes_snps_output_path)
    os.chdir(top_genes_snps_output_path)

    # Output merged data 
    query_output_of.to_csv("top_genes_snps.csv", index=None)


    add_log(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")


'''
######################
###  Bubble Plot #####
######################
'''
def bubbleplot(projectName):

    global projectID 

    # load module 
    ggplot2 = importr('ggplot2')
    dplyr = importr('dplyr',  on_conflict="warn")

    # current path 
    currentPath = os.getcwd()

    # back to root path 
    os.chdir('../../../')
    root_path = os.getcwd()

    # get gwas_lead_snp.txt file 
    input_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME + '/top_genes/'
    if not os.path.exists(input_path): 
        warning = "Please make sure that you have created an input folder with input files"
        add_log(warning)
        add_log('Input path should be: %s' %input_path)
    os.chdir(input_path)

    gwas_lead_snp = r("""
        gwas_lead_snp <- read.csv('sorted_top_genes.csv') %>% 
          select(gene_name, chr, start) %>%
          distinct()
        colnames(gwas_lead_snp) = c('snpsNames', 'chrosome', 'startSites')

        """)
    robjects.globalenv['dataframe'] = gwas_lead_snp


    # annotated metaxcan file path 
    merged_annoated_files_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME + '/annotated_metaxcan_output_files/merged/'
    os.chdir(merged_annoated_files_path)

    data = r("data <- read.csv('merged_annotated_metaxcan_output.csv')")
    robjects.globalenv['dataframe'] = data 

    bubble_plot_output_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME + '/bubble_plot/'
    if not os.path.exists(bubble_plot_output_path): 
         os.makedirs(bubble_plot_output_path)
    os.chdir(bubble_plot_output_path)

    # for loop through each snps site +/- 1000000 bp  
    r("""
        # read dataframe 
        startSites = as.numeric(as.character(gwas_lead_snp$startSites))
        snpsNames = as.character(gwas_lead_snp$snpsNames)
        chrosome = as.numeric(as.character(gwas_lead_snp$chrosome))

        # for loop through each snps site +/- 1000000 bp  
        for (i in 1:length(startSites)) 
        {
            print(paste("BUBBLE PLOTING", ": ", snpsNames[i]))

            # subset data for each snps 
            subData <- subset(data, data$start > startSites[i] - 1000000  & 
                data$start < startSites[i] + 1000000 & data$chr==chrosome[i])
            subData <- subData[order(subData$start),]
            subData <- mutate(subData, z_score=ifelse(subData$zscore > 0, '   +   ', '   -   ')) 
            # Labels = subData$gene_name
            subData$gene_name <- factor(subData$gene_name, levels=subData$gene_name)
            # subData$tissue <- factor(subData$tissue, levels=subData$tissue)
            subData <- subData %>% 
            inner_join(tissue_abbr, by = "tissue")


            # draw plot 
            p <- ggplot(subData, aes(x=subData$gene_name, y=subData$tissue_site_detail_abbr, size=abs(subData$zscore)))
            p + 
            geom_point(aes(colour=z_score)) + 
            scale_color_manual(values=c('blue', 'brown')) +
            scale_size_continuous(guide=FALSE, range=c(0,max(abs(subData$zscore))))+
            # ggtitle(paste('locus: ', snpsNames[i], '(chromosome', chrosome[i], ')')) +
            labs(x='Gene', y='Tissue') + 
            ggtitle(snpsNames[i])+
            # scale_x_discrete(breaks = subData$gene_name, labels=Labels) + 
            theme(axis.text.x = element_text(size=16, face='bold', angle = 90, hjust = 1)) +
            theme(axis.text.y = element_text(size=16, face='bold')) +
            theme(plot.title = element_text(size=18, face='bold')) +
            # theme(axis.title= element_text(size=18, face='bold')) +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
            theme(legend.position = "none") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

            # save plot 
            ggsave(paste(snpsNames[i], '_bubble_plot', '.png', sep=''), width=12, height=12)
        } 
    """)

    add_log(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")


'''
######################
###  Region Plot #####
######################
'''
def regionplot(projectName):

    global projectID 

    # load module 
    ggplot2 = importr('ggplot2')
    dplyr = importr('dplyr',  on_conflict="warn")

    # current path 
    currentPath = os.getcwd()

    # back to root path 
    os.chdir('../../../')
    root_path = os.getcwd()

    # get gwas_lead_snp.txt file 
    input_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME + '/top_genes/'
    if not os.path.exists(input_path): 
        warning = "Please make sure that you have created an input folder with input files"
        add_log(warning)
        add_log('Input path should be: %s' %input_path)
    os.chdir(input_path)

    gwas_lead_snp = r("""
        gwas_lead_snp <- read.csv('sorted_top_genes.csv') %>% 
          select(gene_name, chr, start) %>%
          distinct()
        colnames(gwas_lead_snp) = c('snpsNames', 'chrosome', 'startSites')

        """)
    robjects.globalenv['dataframe'] = gwas_lead_snp

    # annotated metaxcan file path 
    merged_annoated_files_path = root_path + '/output/' + projectName +  "_" + CURRENT_TIME + '/annotated_metaxcan_output_files/merged/'
    os.chdir(merged_annoated_files_path)

    data = r("data <- read.csv('merged_annotated_metaxcan_output.csv')")
    robjects.globalenv['dataframe'] = data 

    bubble_plot_output_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME + '/region_plot/'
    if not os.path.exists(bubble_plot_output_path): 
         os.makedirs(bubble_plot_output_path)
    os.chdir(bubble_plot_output_path)

    r("""
        # read dataframe 
        startSites = as.numeric(as.character(gwas_lead_snp$startSites))
        snpsNames = as.character(gwas_lead_snp$snpsNames)
        chrosome = as.numeric(as.character(gwas_lead_snp$chrosome))

        # loop through snps site +/- 1000,000 bps 
        for (i in 1:length(startSites)) 
        {
        print(paste("REGION PLOTING", ": ", snpsNames[i]))

        # subset data 
        subData <- subset(data, data$start > startSites[i] - 1000000  & data$start < startSites[i] 
            + 1000000 & data$chr==chrosome[i])
        subData$logp <- -log10(subData$pvalue)
        subData <- subData[order(subData$start),]

        subData <- mutate(subData, sig=ifelse(subData$logp > -log10(0.05/nrow(data)), 'Most Sig', 
        ifelse(subData$logp > 5.30103 & subData$logp <= -log10(0.05/nrow(data)), 'Sig', 
        ifelse(subData$logp > -log10(0.05/nrow(subData)) 
            & subData$logp <= 5.30103, 'Less Sig','Not Sig')))) 
        subData$gene_name <- factor(subData$gene_name, levels=subData$gene_name)

        p <- ggplot(subData, aes(x=gene_name, y=logp))
        p + geom_point(aes(colour = sig)) + 
        scale_color_manual(guide=FALSE, values=c('black', 'black', 'black', 'black')) +
            # ggtitle(paste("locus: ", snpsNames[i], '(chromosome', chrosome[i], ')')) + 
        labs(x='Gene', y='-log10(p-value)') +
        ggtitle(snpsNames[i])+
        geom_hline(yintercept = 5.30103, linetype='dashed', color='black') +
        geom_hline(yintercept = -log10(0.05/nrow(data)), color='black') +
        geom_hline(yintercept = -log10(0.05/nrow(subData)), linetype='dashed', color='black') +
        theme(axis.text.x = element_text(size=16, face='bold', angle = 90, hjust = 1)) +
        theme(axis.text.y = element_text(size=16, face='bold')) + 
        theme(plot.title = element_text(size=18, face='bold')) +
        theme(axis.title.x = element_blank(), axis.title.y=element_text(size=18, face='bold')) + 
        theme(legend.position = "none")+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

        # output plot 
        ggsave(paste(snpsNames[i], '_region_plot', '.png',sep=''), width=12, height=12)
        } 
    """) 

    add_log(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")


'''
#########################
###  locuszoom Plot #####
#########################
'''
def locuszoom_plot(projectName):

    dplyr = importr('dplyr',  on_conflict="warn")

    global projectID 

    # current path 
    currentPath = os.getcwd()

    # back to root path 
    os.chdir('../../../')
    root_path = os.getcwd()

    # make batch_locuszoom.txt file 
    input_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME + '/top_genes/'
    if not os.path.exists(input_path): 
        warning = "Please make sure that you have created an input folder with input files"
        add_log(warning)
        add_log('Input path should be: %s' %input_path)
    os.chdir(input_path)

    gwas_lead_snp = r("""
        gwas_lead_snp <- read.csv('sorted_top_genes.csv') %>% 
          select(gene_name, chr, start, end) %>%
          distinct() %>% 
          mutate(flank='1.0MB', run = 'yes', m2zargs= "showAnnot=F")
          colnames(gwas_lead_snp) = c('snp', 'chr', 'start', 'stop', 'flank', 'run', 'm2zargs')
          write.table(gwas_lead_snp, file="batch_locuszoom.txt", sep=" ", quote=FALSE, row.names=F)

        """)


    # write gwas_lead_snp into txt file
    # gwas_lead_snp.to_csv("batch_locuszoom.txt", sep='\t')

    # create a folder to hold locuszoom results 
    locuszoom_plot_path = root_path + '/locuszoom/locuszoom_plots/'
    if not os.path.exists(locuszoom_plot_path): 
        os.makedirs(locuszoom_plot_path)

    # source files from input path 
    # copy the files including plink, run_locuszoom.py and two other .txt from input path into locuszoom program 
    destination = locuszoom_plot_path     # locuszoom plot destination path 
    plink_destination = root_path + '/locuszoom'


    source = os.listdir(input_path)


    for file in source: 
        if file.endswith(".txt"):
            shutil.copy(file, destination)
            msg = "COPYING: the file '%s' into the folder '%s'" % (file, destination) 
            add_log(msg)
    os.system('rm batch_locuszoom.txt')


    # input path 
    input_path = root_path + '/input/'
    if not os.path.exists(input_path): 
        warning = "Please make sure that you have created an input folder with input files"
        add_log(warning)
        add_log('Input path should be: %s' % input_path)
    os.chdir(input_path)

    source = os.listdir(input_path)

    plink_file = 'plink'
    shutil.copy(plink_file, plink_destination)

    msg = "COPYING: the file '%s' into the folder '%s'" % (plink_file, plink_destination)
    add_log(msg)

    for file in source: 
        if file.endswith(".txt"):
            shutil.copy(file, destination)
            msg = "COPYING: the file '%s' into the folder '%s'" % (file, destination) 
            add_log(msg)

    locuszoom_path = root_path + '/MetaXcan-Postprocess/scripts/Python_scripts'
    os.chdir(locuszoom_path)

    locuszoom_scrip_cmd = 'run_locuszoom.py'
    shutil.copy(locuszoom_scrip_cmd, destination)
    msg = "COPYING: the file '%s' into the folder '%s'" % (locuszoom_scrip_cmd, destination)  
    add_log(msg)

    # run locuszoom program  
    os.chdir(destination)
    os.system('./run_locuszoom.py')
    add_log('RUNNING: run_locuszoom.py...')

    # setup locuszoom plot output file path and move all files to this new output folder 
    locuszoom_plot_files_path = root_path + '/output/' + projectName + "_" + CURRENT_TIME + '/looszoom_plot/'
    if not os.path.exists(locuszoom_plot_files_path): 
        os.makedirs(locuszoom_plot_files_path)
    os.chdir(destination)
    
    src = destination
    dst = locuszoom_plot_files_path
    os.system('mv %s %s' % (src, dst))

    msg = "MOVING: all files from the folder '%s' into the folder '%s'" % (src, locuszoom_plot_files_path)
    add_log(msg)


    add_log(datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done!")







