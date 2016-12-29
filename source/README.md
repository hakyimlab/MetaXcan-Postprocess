# Source Scripts

This folder contains a rag tag collection of scripts for plotting and analizing MetaXcan results.

## ProcessMetaXcanresultsFolder

This script will process Metaxcan Results files from a folder and produce QQ plots and Manhattan plots for each of them.

General usage looks like:

```bash
./ProcessMetaXcanResultsFolder.py \
--gencode_file data/gencode.gtf.gz \
 --pattern ".*DGN.*" \
--input_folder data/metaxcan_results \
--output_folder results/images
```

, which will pick any files from `data/metaxcan_results`,
with a `DGN` component in their file name string,
and produce output in `results/images`.

It requires a gencode file such as this:
`ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz`

Please bear in mind hat at the time of this writing, MetaXcan uses a gene annotation that is derived from Gencode release 19, GTEx v6p
```
# GTEx login required:
http://www.gtexportal.org/static/datasets/gtex_analysis_v6/reference/gencode.v19.genes.patched_contigs.gtf.gz
```


# Requirements

- Python 2.7
- scipy, numpy, pandas, rpy2
- R 3 or above
- ggplot, qqman, tidyr, dplyr
