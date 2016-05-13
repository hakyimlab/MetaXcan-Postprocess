# All helper functions can be used through MetaXcan post-processing  

######################
#### Load modules ####
######################

import os, sqlite3 
from datetime import datetime
import uuid as myuuid 

logID = str(myuuid.uuid4())

###########################
#### Global variables  ####
###########################

logfile = None
current_logs = None
LINE =  '-----------------------------------------------------------------------------------------'
SQL_QUERY_PREFIX = "select e.genename, w.rsid from weights w join extra e on w.gene = e.gene where e.genename = '"
SQL_QUERY_PREFIX_DNG = "select gene, rsid from weights where gene = '"
CURRENT_TIME = datetime.now().strftime('%Y-%m-%d-%H.%M.%S') 


###################
#### Log files ####
##################$

# Create a log file path 
def get_log_path(projectName):  

    global logID 

    # current path 
    currentPath = os.getcwd()

    os.chdir('..')
    os.chdir('..')
    os.chdir('..')

    # root path 
    root_path = os.getcwd()

    log_path = root_path + "/log/" + projectName + "_" + CURRENT_TIME + "/"
    if not os.path.exists(log_path): 
        os.makedirs(log_path)
    return log_path

# Create a new log file 
def open_log(projectName):

    global logfile
    global current_logs

    logfile = open(get_log_path(projectName) + datetime.now().strftime('%Y.%m.%d.%H:%M:%S_') + projectName + '.txt', "w")
    current_logs = []

# Add logs 
def add_log(log_str):

    global current_logs

    print (log_str)
    current_logs.append(log_str + '\n')

# Output log 
def write_logs():

    global logfile
    global current_logs

    for index in range (len(current_logs)): 
        logfile.write(current_logs[index])

    current_logs = [] # clear up log 

# pre log 
def pre_message(partName):

    messages = [] 
    messages.append(LINE)
    messages.append('     Getting started MetaXcan postprocessing -- %s' %partName)
    messages.append(LINE)
    for msg in messages:
        add_log(msg)

# Close log 
def finish_log():

    global logfile
    logfile.close()




