
COV_FILE = "forensics/OR4C46_cov.txt"
SNP_FILE = "forensics/OR4C46_snps.txt"
OUTPUT_FILE = "forensics/OR4C46.rds"

cov <- read.csv("forensics/OR4C46_cov.txt")
rownames(cov) <- colnames(cov)

snp <- read.csv("forensics/OR4C46_snps.txt")

x = list(cov=cov, snp=snp)

saveRDS(x,OUTPUT_FILE)