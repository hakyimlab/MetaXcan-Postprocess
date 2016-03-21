# All helper functions can be used through MetaXcan post-processing  

######################
#### Load modules ####
######################
try: 
   import os, sqlite3 
   from datetime import datetime
except Exception as e: 
   hints = "Please download latest python3 and modules at https://www.continuum.io/downloads" 
   add_log(hints)
   add_log(e)

###########################
#### Global variables  ####
###########################

logfile = None
current_logs = None
LINE =  '-----------------------------------------------------------------------------------'
OUTPUT_POSTFIX = '_output_file_metaxcan_postprocessing.csv'
SQL_QUERY_PREFIX = "select e.genename, w.rsid from weights w join extra e on w.gene = e.gene where e.genename = '"
databases = [] 

###################
#### Log files ####
##################$

# Create a log file path 
def get_log_path(filename):  
    log_path = get_current_path() + "/log/" + get_project_name(filename) + "/"
    if not os.path.exists(log_path): 
        os.makedirs(log_path)
    return log_path

# Create a new log file 
def open_log(filename):
    global logfile
    global current_logs

    logfile = open(get_log_path(filename) + get_current_time() + filename, "w")
    current_logs = []

# Add logs 
def add_log(log_str):
    global current_logs

    print (log_str)
    current_logs.append(log_str + '\n')

# Output log 
def finish_metaxcan_postprocessing(filename):
    global logfile
    global current_logs

    for index in range (len(current_logs)): 
        logfile.write(current_logs[index])

    current_logs = [] # clear up log 

# Introduction messages 
def pre_message(partName):
    messages = [] 
    messages.append(LINE)
    messages.append('- Getting started MetaXcan postprocessing -- %s- ' %partName)
    messages.append(LINE)
    for msg in messages:
        # add_log(msg)
        pass 

# Conclusion messages 
def post_message(filename):
    messages = [] 
    messages.append(LINE)
    messages.append('output: %s' % get_out_put_file(filename))
    messages.append(LINE)
    messages.append("Done!\n")
    messages.append('@' + datetime.now().strftime('%H:%M:%S %Y/%m/%d\n'))
    for msg in messages:
        # add_log(msg)
        pass 

# Close log 
def finish_log():
    global logfile

    logfile.close()


####################
#### File paths #### 
####################

# Output file path 
def get_out_path(filename):
    out_path = get_current_path() + "/out/" + get_project_name(filename) + "/"
    if not os.path.exists(out_path): 
        os.makedirs(out_path)
    return out_path

# Output file path with file name 
def get_out_put_file(filename):
    return "%s%s%s"%(get_out_path(filename), filename, OUTPUT_POSTFIX)  

# Input file path with file name 
def get_input_path(filename, tag):  # tag - such as .csv 
    input_path = get_current_path() + "/annotate/"+ filename + tag 
    if not os.path.exists(input_path): 
        warning = "Please make sure that you have an input folder with input files"
        add_log(warning)
    # add_log('Input: %s' %input_path)
    return input_path

# Input file path without file name 
def get_input_path_without_filename():
    input_path = get_current_path() + "/input/"
    if not os.path.exists(input_path): 
        warning = "Please make sure that you have an input folder with input files"
        add_log(warning)
    # add_log('Input: %s' %input_path)
    return input_path

# Database file path 
def get_database_path(filename):
    database_path = get_current_path() + "/databases/"
    if not os.path.exists(database_path): 
        warning = "Please make sure that you have created a databases folder" 
        add_log(warning)
    return database_path 

# Current path 
def get_current_path():
    return os.getcwd() 


######################################
#### Cancer name and current time #### 
######################################

# Current cancer name 
def get_project_name(filename):
    return filename.split("_")[0]

# Current time 
def get_current_time():
    return datetime.now().strftime('%Y.%m.%d.%H:%M:%S_') 


####################
#### Databases #####
####################

# Connect databases 
def connect_database(filename, databases, index):
    databaseName = databases[index]        
    add_log('FETCHING SNPs: ' + databaseName)
    return sqlite3.connect(get_database_path(filename)+databaseName) 


