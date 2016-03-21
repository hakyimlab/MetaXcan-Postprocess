--------------------------------------------------------------------------------

<h1 style="text-align: center;" markdown="1"> MetaXcan Post-processing </h1>
<br> 

## Introduction 
+ Python3 scripts used to postprocess outputs from PrediXcan and MetaXcan

## Prerequisites
+ Compatible with both Python2.7/3, but python3 is recommended.    
+ Download latest Python3 @ <a href="https://www.continuum.io/downloads"> Anaconda </a>. 
+ Open the terminal, and type `alias python='python3.5'` to select python3 if there are multiple python versions installed     

## Installation 
+ Clone this software package, and navigate to the directory (referenced as `<dir>` thereafter) where all package files should be located. These files are described below:  

 + `run.sh` - Helper bash script that launches run_MetaXcanPostprocessing file 
 + `run_MetaXcanPostprocessing.py` - A simple "main()" that loads `MetaXcanPostprocessing.py` and `R2Python.py` file
 + `MetaXcanPostprocessing.py` - This is where all post-processing work will be, such as manhattan plot, qqplot, table of top genes with corresponding snps by tissues  
 + `R2Python.py` - This wraps R packages such as annotables and qqman in python environment using rpy2 
 + `helpers.py` - This contains all helper functions that can be used through post-processing   
 + `__init__.py` - This is a marker file that marks current directory as python package directory 
 + `README.md` - A brief description about this software package 

## Setup 
+ If `run.sh` is not excutable, run ```chmod+x run.sh``` first  to make it executable 
+ Download all database files (*.db) into the directory `<dir>/databases` 
+ Place MetaXcan raw results (*.csv) into this directory `<dir>/input`

## Run 
+ Open the terminal, and execute the `run.sh` script by typing:
 
 ```./run.sh <your project title>``` - your project title should be consise and a single word such as breastcancer, ovarycancer and diabetes. 

## Output 
+ Automatically create two new directories which are used to save log and results: 

 + `<dir>/log/` - Including all logs 

 + `<dir>/out/` - Including pvalue-sorted annotated results with SNPs, QQ-Plots, and Manhattan-Plots for each tissues and all tissues  

 + `<dir>/input/` - Including raw results and annotated results for each tissue 

--------------------------------------------------------------------------------

Copyright (C) 2016 Im lab, Section of Genetics, Department of Medicine, University of Chicago. 
