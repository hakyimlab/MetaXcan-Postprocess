# 
#						 MetaXcan Post-Processing Scripts  
#
#             				        	   Im Lab 
#       		       Section of Genetics, Department of Medicine, University of Chicago 
#          							
#
# 					Part one: Top gene lists
# 					Part two: QQ-plot (tissue vs all tissues)
# 					Part three: QQ-plot (known OC genes vs all genes)
# 					Part four: Manhattan plot (tissue vs all tissues)
# 
# 		
###############################################
#### How to run this script from terminal  ####
###############################################

# install latest python 3 if not installed:
#          https://www.continuum.io/downloads  # this software packages containing many modules we will import, such as pandas, sqlite3 etc. 
#          alias python='python3.5'   # run this in the terminal to select python3.5 if there are multiple python versions installed  
# cd to your working directory containing: 
#          this python script (MetaXcan_postprocessing.py)
#          all database files (*.db)
#          input_file_top_gene_list.csv 
# Generically, it is run as follows in the terminal:
#          python MetaXcan_postprocessing.py 


################################## Part one: Top gene lists with SNPs ####################################

#### A file need to be parpared to run Part one's Script ####

# Prepare "input_file_top_gene_list.csv" file as follow: 
	# column 1 (databases): please keep the order of database names (from top to bottom) with the same order as tissue column names (from column 2 to last column) 
	# column 2 (CrossTissue): gene names 
	# column 3 (Adipose): gene names 
	# .
	# .
	# column n (Breast_mammary): gene names 

		# an example of "input_file_top_gene_list.csv"
		#   | databases  				    | CrossTissue | TW_Adipose-Subcutaneous | TW_Breast-MammaryTissue | TW_Ovary |
		#   | CrossTissue_elasticNet0_0.5.db    	    | CHMP4C 	  | CHMP4C                  | CHMP4C                  | LRRC37A  | 
		#   | TW_Adipose-Subcutaneous_elasticNet0_0.5.db    | CRHR1       | CRHR1                   | CRHR1                   | LRRC37A2 | 
		#   | TW_Breast-MammaryTissue_elasticNet0_0.5.db    | HOXB2       | HOXB2                   | HOXB2                   | PTX3     | 
		#   | TW_Ovary_elasticNet0_0.5.db      		    | HOXB3       | LEKR1                   | HOXB9                   | VEPH1    | 
		#   |					            | HOXD1       | LRRC37A                 | LEKR1                   |          | 
		#   |						    | HOXD3       | LRRC37A2                | LRRC37A                 |          | 
		#   |						    | LEKR1       |                         | LRRC37A2                |          |
		#   |						    | LRRC37A     |                         | WNT3                    |          | 

print ('------------------------------------------------------')
print ('- Getting started, Part one: Top gene list with SNPs -')
print ('------------------------------------------------------')

try: 
   import os, pandas, sqlite3 
except ImportError: 
   raise ImportError ("Please download latest python3 and modules at https://www.continuum.io/downloads ")

def fetchTopGeneList(inputFile="input_file_top_gene_list.csv"):

	#read input file 
	print('Input: input_file_top_gene_list.csv')
	
	try:
		input = pandas.read_csv(inputFile)
	except Exception as e:
		print("Errors in reading:", inputFile)
		print(e)  

	headers = input.dtypes.index 
	tissue_names = headers[1:]

	databases = input[headers[0]]
	databases = databases.dropna()

	print('------------------------------------------------------')
	print("A list of tissues: ")
	gene_lists = [] 
	for i in range(len(tissue_names)):
		print (tissue_names[i])
		genes = input[tissue_names[i]]
		genes = genes.dropna()
		gene_lists.append(genes)
	print('------------------------------------------------------')

	# connect and query through all databases 
	query_output_merged = pandas.DataFrame()
	query_output_list_all = []

	for i in range(len(databases)):
		conn = sqlite3.connect(databases[i])   # Connecting to the database file 
		print("Connect and query database: %s" %databases[i])

		#prefix for query 
		genename_sql_query_prefix = "select e.genename, w.rsid from weights w join extra e on w.gene = e.gene where e.genename = '"

		# get a list of full query 
		full_query_name_list = []
		for k in range(len(gene_lists[i])):
			full_query_name = genename_sql_query_prefix + gene_lists[i][k] + "'"
			full_query_name_list.append(full_query_name)

		# query database by looping through all genes
		query_output_list = []
		for m in range(len(full_query_name_list)):
			query_output = pandas.read_sql(full_query_name_list[m], conn, index_col=None)
			query_output['Tissue'] = tissue_names[i] #add correspinding tissue names to the new 'Tissue' column  
			query_output_list.append(query_output) 

		query_output_of = pandas.concat(query_output_list, axis = 0) #merge output data from current database 

		query_output_list_all.append(query_output_of)
		
		conn.close()    # Closing the connection to the database file

	query_output_merged = pandas.concat(query_output_list_all, axis=0)  #merge all output data from all databases 

	#output data 
	query_output_merged.to_csv("output_file_top_gene_list.csv", index=None)

	print('------------------------------------------------------')
	print('output: output_file_top_gene_list.csv')
	print('------------------------------------------------------')
	print("Done!\n")


if __name__ == "__main__":
	fetchTopGeneList()







