# This is where all post-processing work will be, 
# such as manhattan plot, qqplot, table of top genes with corresponding snps by tissues. 

try: 
    from helpers import * 
    from datetime import datetime 
    import os, pandas, sqlite3, glob
except Exception as e: 
    warning = "Please download latest python3 and modules at https://www.continuum.io/downloads "
    print(warning)
    print(e)
    add_log(warning)
    add_log(e)

############################################
#### Part one: Top gene lists with SNPs ####
############################################

class TopGeneListSnps(object):

    # Constructor 
    def __init__(self, filename):
        # Declare all instance variables 
        self.tissue_names = []
        self.databases = []
        self.gene_lists = [] 

        self.query_output_list_all = []
        self.query_output_merged = pandas.DataFrame()

        # Get all information about databases, tissues, and corresponding genes 
        self.readFile(filename)  

    # Read input file 
    def readFile(self, filename):
        # Log introduction messages 
        pre_message("Part one: Top gene lists with snp ")

        # Set up input file path  
        input_file_path = get_input_path(filename, ".csv")

        # Read input file 
        try:
            input = pandas.read_csv(input_file_path)
            if not os.path.exists(input_file_path):
                msg = 'please double check there is input file in the path: %s' %input_file_path 
                print(msg)
                add_log(msg)
        except Exception as e:
            msg = "Errors in reading: %s" %input_file_path 
            print(msg)
            add_log(msg)
            print(e) 
            add_log(e) 

        # Get tissue names (2nd to last)
        headers = input.dtypes.index 
        self.tissue_names = headers[1:]

        # Get database names 
        self.databases = input[headers[0]]
        self.databases = self.databases.dropna() # Remove empty elements  

        print(LINE)  
        add_log(LINE)
        msg_tissue = "A list of tissues: "
        print(msg_tissue)
        add_log(msg_tissue)


        # Fetch a list of genes from input file  
        for i in range(len(self.tissue_names)):

            print (self.tissue_names[i])  # Print tissue names 
            add_log(self.tissue_names[i])

            genes = input[self.tissue_names[i]]
            genes = genes.dropna()    # Remove empty elements
            self.gene_lists.append(genes)

        print(LINE)
        add_log(LINE)  

    # Fetch top gene list with snps from databases 
    def fetchTopGeneList(self, filename):
        # Loop through databases 
        for i in range(len(self.databases)):
            # Connect databases 
            conn = connect_database(filename, self.databases, i)

            # Get a list of full query 
            full_query_name_list = []
            for k in range(len(self.gene_lists[i])):
                full_query_name = SQL_QUERY_PREFIX + self.gene_lists[i][k] + "'"
                full_query_name_list.append(full_query_name)

            # Looping through all genes
            query_output_list = []
            for m in range(len(full_query_name_list)):
                query_output = pandas.read_sql(full_query_name_list[m], conn, index_col=None)

                # Add correspinding tissue names to the new 'Tissue' column 
                query_output['Tissue'] = self.tissue_names[i] 
                query_output_list.append(query_output) 

            # Merge output data 
            query_output_of = pandas.concat(query_output_list, axis = 0) 
            self.query_output_list_all.append(query_output_of)
            
            # Close database
            conn.close()    

    # Output results 
    def outputTopGeneList(self, filename): 
        # Merge all output data
        self.query_output_merged = pandas.concat(self.query_output_list_all, axis=0) 

        # Output merged data 
        self.query_output_merged.to_csv(get_out_put_file(filename), index=None)

        # Log conclusion messages  
        post_message(filename)

    # Output logs 
    def logTopGeneList(self, filename):
        finish_metaxcan_postprocessing(filename) 





