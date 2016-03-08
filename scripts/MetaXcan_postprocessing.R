# 
#						           MetaXcan Post-Processing Scripts  
#
#             				        	   Im Lab 
#              Section of Genetics, Department of Medicine, University of Chicago 
#          							
#
# 						Part one: Top gene lists with SNPs 
# 						Part two: QQ-plot (tissue vs all tissues)
# 						Part three: QQ-plot (known OC genes vs all genes)
# 						Part four: Manhattan plot (tissue vs all tissues)
# 		

####################################################
#### How to run this script from R command line ####
####################################################
# install R if not installed 
#          https://cran.cnr.berkeley.edu
# cd to your working directory containing: 
#          this R script (MetaXcan_postprocessing.R)
#          all database files (*.db)
#          input_file_top_gene_list.csv 
# open R in the terminal by typing 'R' and then hit 'return' key 
# Generically, it is run as follows in the R command line: 
#           source('MetaXcan_postprocessing.R')


################################## Part one: Top gene lists with SNPs ####################################

##################################################################
####    A file need to be parpared to run Part one's Script   ####
##################################################################
# Prepare "input_file_top_gene_list.csv" file as follow: 
	# column 1 (databases): please keep the order of database names (from top to bottom) with the same order as tissue column names (from column 2 to last column) 
	# column 2 (CrossTissue): gene names 
	# column 3 (Adipose): gene names 
	# .
	# .
	# column n (Breast_mammary): gene names 

		# an example of "input_file.csv"
		#	| databases  									| CrossTissue | TW_Adipose-Subcutaneous | TW_Breast-MammaryTissue | TW_Ovary |
		#   | CrossTissue_elasticNet0_0.5.db    			| CHMP4C 	  |	CHMP4C                  | CHMP4C                  | LRRC37A  | 
		#	| TW_Adipose-Subcutaneous_elasticNet0_0.5.db    | CRHR1       | CRHR1                   | CRHR1                   | LRRC37A2 | 
		#	| TW_Breast-MammaryTissue_elasticNet0_0.5.db    | HOXB2       | HOXB2                   | HOXB2                   | PTX3     | 
		#	| TW_Ovary_elasticNet0_0.5.db      				| HOXB3       | LEKR1                   | HOXB9                   | VEPH1    | 
		#   |												| HOXD1       | LRRC37A                 | LEKR1                   |          | 
		#	|												| HOXD3       | LRRC37A2                | LRRC37A                 |          | 
		#	|												| LEKR1       |                         | LRRC37A2                |          |
		#	|												| LRRC37A     |                         | WNT3                    |          | 


############################
#### 1. Load libraries  ####
############################
message ('------------------------------------------------------')
message ('- Getting started, Part one: Top gene list with SNPs -')
message ('------------------------------------------------------')

#install and load RSQLite library 
if(!require("RSQLite")) {
  message("Please install the RSQLite plackage")
  install.packages("RSQLite")
}
require(RSQLite)

#############################
#### 2. Read input file  ####
#############################
# read input file.csv 
filename = "input_file_top_gene_list.csv"  
input <- read.csv(filename)
message('Input:', filename)

# get tissue names (2nd to last)
column_names <- colnames(input)
tissue_names <- column_names[-1]

# get database name
databases <- as.vector(input[1])
databases <- databases[databases != ""] #remove empty elements 

# get each gene list 
message('------------------------------------------------------')
message("A list of tissues: ")
gene_lists <- list()
for (i in 1:length(tissue_names))
{
	 message(tissue_names[i])
	 genes <- as.vector(input[i+1])
	 genes <- genes[genes != ""] #remove empty elements 
	 gene_lists[[i]] <- genes
}
message('------------------------------------------------------')


##############################################
#### 3. Connect sqlite and query database ####
##############################################
sqlite <- dbDriver("SQLite")
query_output_merged <- data.frame()

#loop through all databases 
for (i in 1: length(databases))
{
	current_databases <- databases[i] 
	db <- dbConnect(sqlite, current_databases) # open database connection
	message ('Connect and query database:', current_databases)

	#prefix for query 
	genename_sql_query_prefix <- "select e.genename, w.rsid from weights w join extra e on w.gene = e.gene where e.genename = '"

	#length of genes_list 
	length_gene_list <- length(gene_lists[[i]]); 

	# get a list of full query 
	full_query_name_vector <- 0 
	for (k in 1:length_gene_list) {
		full_query_name <- paste(genename_sql_query_prefix, gene_lists[[i]][k], c("'"), sep="")
		full_query_name_vector[k] <- full_query_name
	}

	# query database by looping through all genes
	for (m in 1:length(full_query_name_vector))
	{
		query_output <- dbGetQuery(db, full_query_name_vector[m])
		query_output <- transform(query_output, Tissue = tissue_names[i]) #add correspinding tissue names to the new 'Tissue' column 
	    # print(query_output)
		query_output_merged <- rbind(query_output_merged, query_output) #merge output data 
	}
	
	dbDisconnect(db) # close current database connection 
}

########################
#### 4. Output data ####
########################
write.csv(query_output_merged, file="output_file_top_gene_list.csv", row.names=F)
message('------------------------------------------------------')
message('output: output_file_top_gene_list.csv')
message('------------------------------------------------------')
message("Done!\n")




