--------------------------------------------------------------------------------

<h1 style="text-align: center;" markdown="1"> MetaXcan Post-processing </h1>

## Introduction 
+ A pipeline used to postprocess MetaXcan results 

## Prerequisites
+  <a rel="nofollow" class="external text" href="http://www.python.org/download/">Python 2.7+</a> (please do <b> not </b> download the 3.0 branch!)
+  <a rel="nofollow" class="external text" href="http://www.r-project.org/">R 3.0+</a>
+  <a rel="nofollow" class="external text" href="http://rpy2.readthedocs.io/en/version_2.7.x/"> rpy2 </a>
+  <a rel="nofollow" class="external text" href="https://github.com/stephenturner/annotables#how"> annotables </a>
+  <a rel="nofollow" class="external text" href="https://github.com/hadley/dplyr"> dplyr </a>
+  <a rel="nofollow" class="external text" href="https://github.com/stephenturner/qqman"> qqman </a>
+  <a rel="nofollow" class="external text" href="https://github.com/hadley/ggplot2"> ggplot2 </a>

## Installation and Setup 
 + Download this pipline as `<dir>/MetaXcan-Postprocess`:
 + Navigate to the `src`.  
 + `MetaXcanPostprocessing.py` - This is where all post-processing computation will be. 
 + `run_locuszoom.py` - A script that runs standalone software `locuszoom`  
 + `helpers.py` - This contains all helper functions   
 + `__init__.py` - This is a marker file that marks current directory as python package directory 
 + `README.md` - A brief description about this pipline 
 + <a href = "http://genome.sph.umich.edu/wiki/LocusZoom_Standalone"> Download LocusZoom (37.5GB) </a> as `<dir>/locuszoom`
 + Create a new fold as `<dir>/input`, and include the following files:  
  + <a href = "http://hakyimlab.org/predictdb/"> Prediction models </a> (`*.db`)
  + <a href = "http://pngu.mgh.harvard.edu/~purcell/plink/"> plink </a> (`plink`)
  + MetaXcan outputs (`*.csv`) 
  + Input file for locuszoom plot - `gwas_snp.txt`. The text file should be tab-delimited and prepared as described in <a href = "http://genome.sph.umich.edu/wiki/LocusZoom_Standalone"> here </a>. For example: 
  
      MarkerName |	P-value
      ---- | -----
      rs1	|  0.423
      rs2 |	1.23e-04
      rs3 |	9.4e-390

## Run 
+ Open the terminal, and execute:
 ```python MetaXcanPostprocessing.py <your project title>``` 

## Output 
+ two new folders which are used to save log information and results: 
 + `<dir>/log/` - all log files 
 + `<dir>/out/` - the following tables and figures 
    + annotation `.csv`
    + manhattan plot `.png`
    + qqplot `.png`
    + region plot `.png`
    + bubble plot `.png`
    + locuszoom plot `.pdf`
    + top genes `.csv`
    + top genes with snps `.csv` 


--------------------------------------------------------------------------------