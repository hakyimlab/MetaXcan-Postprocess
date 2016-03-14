--------------------------------------------------------------------------------

<h1 style="text-align: center;" markdown="1"> MetaXcan Post-processing </h1>
<br> 

## Introduction 
+ Python3 scripts used to postprocess outputs from PrediXcan and MetaXcan

## Prerequisites  
+ Download latest Python3 @ <a href="https://www.continuum.io/downloads"> Anaconda </a>. 
+ Open the terminal, and type `alias python='python3.5'` to select python3 if there are multiple python versions installed     

## Installation 
+ Clone this software package, and navigate to the directory (referenced as `<dir>` thereafter) where all package files should be located. These files are described below:  

 + `run.sh` - Helper bash script that launches run_MetaXcanPostprocessing file 
 + `run_MetaXcanPostprocessing.py` - A simple "main()" that loads MetaXcanPostprocessing file
 + `MetaXcanPostprocessing.py` - This is where all post-processing work will be, such as manhattan plot, qqplot, table of top genes with corresponding snps by tissues  
 + `helpers.py` - This contains all helper functions that can be used through post-processing   
 + `__init__.py` - This is a marker file that marks current directory as python package directory 
 + `README.md` - A brief description about this software package 



## Setup 
+ If `run.sh` is not excutable, run ```chmod+x run.sh``` first  to make it executable 
+ Create a new directory `<dir>/databases`, download all database files (*.db) into this directory  
+ Create a new directory `<dir>/input`, place all data input files into this directory  

## Input 

+ Prepare data input files as follow: 
	+ column 1 (databases): please keep the order of database names (from top to bottom) with the same order as tissue column names (from column 2 to last column) 
	+ column 2 (CrossTissue): gene names 
	+ column 3 (Adipose): gene names 
	+ .
	+ .
	+ column n (Breast_mammary): gene names 
+ For example, **ovarycancer\_input\_file\_metaxcan\_postprocessing.csv** (The first word in the filename should be the name of cancer type without any space, e.g. breastcancer, lungcancer etc) 

|  | databases | CrossTissue | TW_Adipose-Subcutaneous | TW_Breast-MammaryTissue | TW_Ovary |
|----|---------|-------|-----------------|---------------|-----------|
| 1 | CrossTissue\_elasticNet0\_0.5.db  | CHMP4C | CHMP4C  | CHMP4C  | LRRC37A |
| 2 | TW\_Adipose\-Subcutaneous\_elasticNet0\_0.5.db | CRHR1  | CRHR1  | CRHR1  | LRRC37A2 |
| 3 | TW\_Breast\-MammaryTissue\_elasticNet0\_0.5.db |  HOXB2  |  HOXB2  |  HOXB9  | PTX3 |
| 4 | TW\_Ovary\_elasticNet0\_0.5.db  |HOXB3 |LRRC37A | LRRC37A | VEPH1|
| 5 | | HOXD1 |LRRC37A2  |HOXB2   |  |
| 6 |  |HOXD3  |  | LRRC37A2 |  |
| 7 |   | LEKR1 |  |  WNT3  |  |
| 8 |  | LRRC37A | | LEKR1 |  |
| 9 |  | |  | VEPH1  |  |
| 10 |  |  |  | PTX3 |  |

## Run 
+ Open the terminal, and execute the `run.sh` script by typing:
 
 ```./run.sh <your input file name>``` 

## Output 
+ Automatically create two new directories: 

 + `<dir>/log/` - For log

 + `<dir>/out/` - For results.

--------------------------------------------------------------------------------

Copyright (C) 2016 Im lab, Section of Genetics, Department of Medicine, University of Chicago. 
