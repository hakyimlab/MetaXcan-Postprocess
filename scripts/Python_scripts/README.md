--------------------------------------------------------------------------------

<h1 style="text-align: center;" markdown="1"> MetaXcan Post-processing </h1>

| Version | Description     |
| --------|:---------------:|
| V1.0    | Initial version |
| V1.1    | Add bubble and region plot |
| V1.2    | Integrate LocusZoom plot |

## Introduction 
+ Python3 scripts used to postprocess outputs from PrediXcan and MetaXcan

## Prerequisites
+ Compatible with both Python2.7/3, but python3 is recommended.    
+ Download latest Python3 @ <a href="https://www.continuum.io/downloads"> Anaconda </a>. 
+ Open the terminal, and type `alias python='python3.x'` to select python3 if there are multiple python versions installed     

## Installation and Setup 
+ Navigate to the directory (referenced as `<dir>` thereafter) where all package files should be located. 
+ Download this pipline in `<dir>/MetaXcan-Postprocess`:  
 + `run.sh` - Helper bash script that launches run_MetaXcanPostprocessing file 
 + `run_MetaXcanPostprocessing.py` - A simple "main()" that loads `MetaXcanPostprocessing.py` and `run_locuszoom` file
 + `MetaXcanPostprocessing.py` - This is where all post-processing work will be, such as output annotation, manhattan plot, qqplot, retion plot, bubble plot, locuszoom plot, top genes with and without corresponding SNPs.   
 + `helpers.py` - This contains all helper functions that can be used through post-processing   
 + `__init__.py` - This is a marker file that marks current directory as python package directory 
 + `README.md` - A brief description about this pipline 
+ Download standalone LocusZoom <a href = "http://genome.sph.umich.edu/wiki/LocusZoom_Standalone"> here </a> in `<dir>/locuszoom`
+ Create a new fold as `<dir>/input`, and add the following files (no subfold, only individual files). 
  + Download prediction models (`*.db`) <a href = "http://hakyimlab.org/predictdb/"> Prediction Models </a>
  + Download plink (`plink`, discard the rest) @ <a href = "http://pngu.mgh.harvard.edu/~purcell/plink/"> plink </a> 
  + Add outputs from PrediXcan or MetaXcan analysis (`*.csv`) 
  + Add locuszoom input file (`gwas_snp.txt`, don't use differen file name). The text file (tab-delimited) should be prepared as described in <a href = "http://genome.sph.umich.edu/wiki/LocusZoom_Standalone"> here </a>
  + Add batch-model file (`batch_locuszoom.txt`, don't use differen file name). The text file (tab-delimited) shold be prepared as described in <a href = "http://genome.sph.umich.edu/wiki/LocusZoom_Standalone"> here </a> 


## Run 
+ Open the terminal, and execute the `run.sh` script by typing:
 ```./run.sh <your project title>``` - your project name should be consise and a single word such as breastcancer, ovarycancer and diabetes. 

## Output 
+ Automatically create two new directories which are used to save log and results: 
 + `<dir>/log/` - Including all log files 
 + `<dir>/out/` - Including output annotation `.csv` files, manhattan plot `.png` files, qqplot `.png` files, region plot `.png` files, bubble plot `.png` files, locuszoom plot `.pdf` files, top genes with and without corresponding SNPs `.csv`files


--------------------------------------------------------------------------------

2016 Im lab, Section of Genetics, Department of Medicine, University of Chicago. 
