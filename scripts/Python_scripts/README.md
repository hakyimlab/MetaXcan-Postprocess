--------------------------------------------------------------------------------

<h1 style="text-align: center;" markdown="1"> MetaXcan Post-processing </h1>
<br> 

| Version | Description     |
| --------|:---------------:|
| V1      | Initial version |
| V2      | Integrate LocusZoom |

## Introduction 
+ Python3 scripts used to postprocess outputs from PrediXcan and MetaXcan

## Prerequisites
+ Compatible with both Python2.7/3, but python3 is recommended.    
+ Download latest Python3 @ <a href="https://www.continuum.io/downloads"> Anaconda </a>. 
+ Open the terminal, and type `alias python='python3.x'` to select python3 if there are multiple python versions installed     

## Installation 
+ Navigate to the directory (referenced as `<dir>` thereafter) where all package files should be located. 
+ Download this pipline as `<dir>/MetaXcan-Postprocess`:  
 + `run.sh` - Helper bash script that launches run_MetaXcanPostprocessing file 
 + `run_MetaXcanPostprocessing.py` - A simple "main()" that loads `MetaXcanPostprocessing.py` and `run_locuszoom` file
 + `MetaXcanPostprocessing.py` - This is where all post-processing work will be, such as output annotation, manhattan plot, qqplot, retion plot, bubble plot, locuszoom plot, top genes with and without corresponding SNPs.   
 + `helpers.py` - This contains all helper functions that can be used through post-processing   
 + `__init__.py` - This is a marker file that marks current directory as python package directory 
 + `README.md` - A brief description about this pipline 
+ Download standalone LocusZoom @ <a href = "http://genome.sph.umich.edu/wiki/LocusZoom_Standalone"> Locuszoom Plot </a> as `<dir>/locuszoom`
+ Create a new fold as `<dir>/input`. 
  + Download prediction models @ <a href = "http://hakyimlab.org/predictdb/"> Prediction Models </a>
  + Download plink (only plink, discard the rest) @ <a href = "http://pngu.mgh.harvard.edu/~purcell/plink/"> plink </a> 
  + Add PrediXcan or MetaXcan outputs 
  + Add locuszoom input file (gwas_snp.txt). The text file (Tab-delimited) should be prepared as follow: 
	| MarkerName | P-value |
	| --- | --- |
	| rs1 | 0.983|
	| rs1 | 1.83e-09 |
	| rs1 | 2.44e-08 |
  + Add batch-model file (batch_locuszoom.txt). The text file (Tab-delimited) shold be prepared as the following format: 
    | snp | chr | start | stop | flank | run | m2zargs |
	| --- | --- | --- | --- | --- | --- | --- |
	| rs7983146 | 2 | 1208977889 | 1298977889 | 1.25MB | yes | title = "snp rs7983146" | 
	| STAT5 | 7 | 1408977889 | 1498977889 | 800kb | yes | title = "STAT5" | 
	| rs9983148 | NA | NA | NA | 2.25MB | yes | title = "unknow snp" | 

## Setup 
+ If `run.sh` is not excutable, run ```chmod+x run.sh``` first  to make it executable
+ The `<dir>/databases` folder has one demo database file, and the `<dir>/input` fold has one demo corresponding MetaXcan raw results. You can run demo with these two files. 
+ Download <a href = "https://app.box.com/s/gujt4m6njqjfqqc9tu0oqgtjvtz9860w"> database files </a> into the directory `<dir>/databases` 
+ Move all your MetaXcan raw results (*.csv) into this directory `<dir>/input`

## Run 
+ Open the terminal, and execute the `run.sh` script by typing:
 
 ```./run.sh <your project title>``` - your project title should be consise and a single word such as breastcancer, ovarycancer and diabetes. 

## Output 
+ Automatically create two new directories which are used to save log and results: 

 + `<dir>/log/` - Including all logs 

 + `<dir>/out/` - Including output annotation `.csv` file, manhattan plot `.png`, qqplot `.png`, region plot `.png`, bubble plot `.png`, locuszoom plot `.pdf`, top genes with and without corresponding SNPs `.csv`


--------------------------------------------------------------------------------

2016 Im lab, Section of Genetics, Department of Medicine, University of Chicago. 
