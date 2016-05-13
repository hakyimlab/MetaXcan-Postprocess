--------------------------------------------------------------------------------

<h1 style="text-align: center;" markdown="1"> MetaXcan Post-processing </h1>

## Introduction 
+ Python3 scripts used to postprocess outputs from PrediXcan and MetaXcan

## Prerequisites
+ Compatible with both Python2.7/3, but python3 is recommended.    
+ Download latest Python3 version @ <a href="https://www.continuum.io/downloads"> Anaconda </a>.  
+ Open the terminal, and type `alias python='python3.5'` to select python3 if there are multiple python versions installed 
+ Install the following python or R libraries: 
   + rpy2 python module: `pip install rpy2`. More detail @ <a href="http://rpy2.readthedocs.io/en/version_2.7.x/"> rpy2 </a>
   + annotables R package: `install.packages("devtools")`, and `devtools::install_github("stephenturner/annotables")`. More detail @ <a href="https://github.com/stephenturner/annotables#how"> annotables </a>
   + dplyr R package: `install.packages('dplyr')`. More detail @ <a href="https://github.com/hadley/dplyr"> dplyr </a>
   + qqman R package: 'install.packages("qqman")`. More detail @ <a href="https://github.com/stephenturner/qqman"> qqman </a>

## Installation and Setup 
+ Navigate to the directory (referenced as `<dir>` thereafter) where all package files should be located. 
+ Download this pipline as `<dir>/MetaXcan-Postprocess`:  
 + `run.sh` - Helper bash script that launches run_MetaXcanPostprocessing file 
 + `run_MetaXcanPostprocessing.py` - A simple "main()" that loads `MetaXcanPostprocessing.py` and `run_locuszoom` file
 + `MetaXcanPostprocessing.py` - This is where all post-processing work will be, such as output annotation, manhattan plot, qqplot, retion plot, bubble plot, locuszoom plot, top genes with and without corresponding SNPs.   
 + `helpers.py` - This contains all helper functions that can be used through post-processing   
 + `__init__.py` - This is a marker file that marks current directory as python package directory 
 + `README.md` - A brief description about this pipline 
+ Download standalone LocusZoom <a href = "http://genome.sph.umich.edu/wiki/LocusZoom_Standalone"> here </a> as `<dir>/locuszoom`
+ Create a new fold as `<dir>/input`, and add the following files (don't creat subfolder, only individual files). 
  + Download prediction models (`*.db`) <a href = "http://hakyimlab.org/predictdb/"> Prediction Models </a>
  + Download plink (only `plink` file, discard the rest) @ <a href = "http://pngu.mgh.harvard.edu/~purcell/plink/"> plink </a> 
  + Add outputs from PrediXcan or MetaXcan analysis (`*.csv`) 
  + Add locuszoom input file (`gwas_snp.txt`, other names won't work). The text file (tab-delimited) should be prepared as described in <a href = "http://genome.sph.umich.edu/wiki/LocusZoom_Standalone"> here </a>
      MarkerName |	P-value
      ---- | -----
      rs1	|  0.423
      rs2 |	1.23e-04
      rs3 |	9.4e-390

  + Add batch-model file (`batch_locuszoom.txt`, other names won't work). The text file (tab-delimited) shold be prepared as described in <a href = "http://genome.sph.umich.edu/wiki/LocusZoom_Standalone"> here </a> 
    snp | chr | start | stop | flank | run | m2zargs
    ---- | --- | ---- | ---- | ----- | ---- | ------
    rs7983146 | NA | NA | NA | 500kb | yes | title="My favorite SNP"
    TCF7L2 | NA | NA | NA | 1.25MB | yes | title="TCF7L2 Region" showRecomb=F
    rs7957197 | 12 | 119503590 | 120322280 | NA | yes | showAnnot=F


## Run 
+ Open the terminal, and execute the `run.sh` script by typing:
 ```./run.sh <your project title>``` - your project name should be consise and a single word e.g. breastcancer, ovarycancer, diabetes etc

## Output 
+ Automatically create two new directories which are used to save log and results: 
 + `<dir>/log/` - Including all log files 
 + `<dir>/out/` - Including output annotation `.csv` files, manhattan plot `.png` files, qqplot `.png` files, region plot `.png` files, bubble plot `.png` files, locuszoom plot `.pdf` files, top genes with and without corresponding SNPs `.csv`files


--------------------------------------------------------------------------------

2016 Im lab, Section of Genetics, Department of Medicine, University of Chicago. 
