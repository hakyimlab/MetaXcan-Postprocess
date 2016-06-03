--------------------------------------------------------------------------------

<h1 style="text-align: center;" markdown="1"> MetaXcan Post-processing </h1>

## Introduction 
+ A pipeline used to postprocess outputs from PrediXcan and MetaXcan

## Prerequisites
+  <a rel="nofollow" class="external text" href="http://www.python.org/download/">Python 2.7+</a> (do <b> not </b> download the 3.0 branch!)
+  <a rel="nofollow" class="external text" href="http://www.r-project.org/">R 3.0+</a>
+  <a rel="nofollow" class="external text" href="http://rpy2.readthedocs.io/en/version_2.7.x/"> rpy2 </a>
+  <a rel="nofollow" class="external text" href="http://rpy2.readthedocs.io/en/version_2.7.x/"> annotables </a>
+  <a rel="nofollow" class="external text" href="http://rpy2.readthedocs.io/en/version_2.7.x/"> dplyr </a>
+  <a rel="nofollow" class="external text" href="http://rpy2.readthedocs.io/en/version_2.7.x/"> qqman </a>
+  <a rel="nofollow" class="external text" href="http://rpy2.readthedocs.io/en/version_2.7.x/"> ggplot2 </a>

## Installation and Setup 
+ Navigate to the directory (referenced as `<dir>` thereafter) where all package files should be located. 
+ Download this pipline as `<dir>/MetaXcan-Postprocess`:  
 + `run.sh` - Helper bash script that launches run_MetaXcanPostprocessing file 
 + `run_MetaXcanPostprocessing.py` - A simple "main()" that loads `MetaXcanPostprocessing.py` and `run_locuszoom` file
 + `MetaXcanPostprocessing.py` - This is where all post-processing work will be, such as output annotation, manhattan plot, qqplot, retion plot, bubble plot, locuszoom plot, top genes with and without corresponding SNPs. 
 + `run_locuszoom.py` - A script that runs standalone software `locuszoom`  
 + `helpers.py` - This contains all helper functions that can be used through post-processing   
 + `__init__.py` - This is a marker file that marks current directory as python package directory 
 + `README.md` - A brief description about this pipline 
+ <a href = "http://genome.sph.umich.edu/wiki/LocusZoom_Standalone"> Download LocusZoom (37.5GB) </a> as `<dir>/locuszoom`
+ Create a new fold as `<dir>/input`, and add the following files. 
  + prediction models (`*.db`) <a href = "http://hakyimlab.org/predictdb/"> Prediction Models </a>
  + plink (only `plink` file) @ <a href = "http://pngu.mgh.harvard.edu/~purcell/plink/"> plink </a> 
  + PrediXcan or MetaXcan outputs (`*.csv`) 
  + Input files for locuszoom plot (`gwas_snp.txt`). The text file (tab-delimited, and use the exact same file and column names) should be prepared as described in <a href = "http://genome.sph.umich.edu/wiki/LocusZoom_Standalone"> here </a>. For example: 
  
      MarkerName |	P-value
      ---- | -----
      rs1	|  0.423
      rs2 |	1.23e-04
      rs3 |	9.4e-390

## Run 
+ Open the terminal, and execute the `run.sh` (`<dir>/MetaXcan-Postprocess/scripts/Python_scripts/`) by typing:
 ```./run.sh <your project title>``` - `your project title` should be consise and a single word e.g. breastcancer, ovarycancer, diabetes etc

## Output 
+ two new folders which are used to save log information and results: 
 + `<dir>/log/` - Including all log files 
 + `<dir>/out/` - Including the following tables and figures 
    + annotation `.csv`
    + manhattan plot `.png`
    + qqplot `.png`
    + region plot `.png`
    + bubble plot `.png`
    + locuszoom plot `.pdf`
    + top genes `.csv`
    + top genes with snps `.csv` 


--------------------------------------------------------------------------------

2016 Im lab, Section of Genetics, Department of Medicine, University of Chicago. 
