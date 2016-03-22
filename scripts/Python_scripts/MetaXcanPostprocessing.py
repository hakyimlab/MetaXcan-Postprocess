# This is where all post-processing work will be, 
# such as manhattan plot, qqplot, table of top genes with corresponding snps by tissues. 

try: 
    from helpers import * 
    from datetime import datetime 
    import os, pandas, sqlite3, glob
except Exception as e: 
    warning = "Please download latest python3 and appropriate modules at https://www.continuum.io/downloads "
    add_log(warning)
    add_log(e)

############################################
#### Part one: Top gene lists with SNPs ####
############################################

class TopGeneListSnps(object):

    # Constructor 
    def __init__(self, filename, fullFilename):
        # Declare all instance variables 
        self.databases = []
        self.gene_lists = [] 

        self.query_output_list_all = []
        self.query_output_merged = pandas.DataFrame()
        self.input = pandas.DataFrame() 


        # Get all information about databases, tissues, and corresponding genes 
        self.readFile(filename)  

    # Read input file 
    def readFile(self, filename):
        # Log introduction messages 
        pre_message("Part two: Top gene lists with snp ")
        currentPath = get_current_path()

        # Set up current and input file path  
        input_file_path = currentPath + '/input/' + 'annotate/' + filename + '.csv'
        # print(input_file_path)

        # Read input file and fetch gene list 
        try:
            self.input = pandas.read_csv(input_file_path)
            if not os.path.exists(input_file_path):
                msg = 'please double check there is input file in the path: %s' %input_file_path 
                add_log(msg)
            self.gene_lists = self.input['gene_name']  
            # print(self.gene_lists)
        except Exception as e:
            msg = "Errors in reading: %s" %input_file_path 
            add_log(msg)
            add_log(e) 

        # Get database names 
        dataPath = get_database_path(filename) 
        os.chdir(dataPath)
        dbFileList = glob.glob("*.db")

        for dbFilename in dbFileList:
           self.databases.append(dbFilename)

        os.chdir(currentPath)    

    # Fetch top gene list with snps from databases 
    def fetchTopGeneList(self, filename):
        # Loop through databases 
        for i in range(len(self.databases)):
            if filename[:-10] in self.databases[i]: 
                # Connect databases 
                conn = connect_database(filename[:-10], self.databases, i)

                # Get a list of full query 
                full_query_name_list = []
                for k in range(len(self.gene_lists)):
                    # print(self.databases[i])
                    if self.databases[i] == 'DGN-WB-unscaled_0.5.db':
                        full_query_name = SQL_QUERY_PREFIX_DNG + self.gene_lists[k] + "'"
                    else: 
                        full_query_name = SQL_QUERY_PREFIX + self.gene_lists[k] + "'"
                    full_query_name_list.append(full_query_name)

                # Looping through all genes
                query_output_list = []
                for m in range(len(full_query_name_list)):
                    query_output = pandas.read_sql(full_query_name_list[m], conn, index_col=None)

                    # Add correspinding parameters to the new output file 
                    query_output['tissue'] = filename[:-10]
                    query_output['pvalue'] =  self.input['pvalue'][m]
                    query_output['zscore'] =  self.input['zscore'][m]
                    query_output['model_n'] = self.input['model_n'][m]
                    query_output['chr'] = self.input['chr'][m]  
                    query_output['start'] = self.input['start'][m] 
                    query_output_list.append(query_output) 

                # Merge output data 
                query_output_of = pandas.concat(query_output_list, axis = 0) 
                self.query_output_list_all.append(query_output_of)
                
                # Close database
                conn.close()    

    # Output results 
    def outputTopGeneList(self, fullFilename): 
        # Merge all output data
        self.query_output_merged = pandas.concat(self.query_output_list_all, axis=0) 

        # Output merged data 
        self.query_output_merged.to_csv(get_out_put_file(fullFilename), index=None)

        # Log conclusion messages  
        post_message(fullFilename)

    # Output logs 
    def logTopGeneList(self, filename):
        finish_metaxcan_postprocessing(filename) 





