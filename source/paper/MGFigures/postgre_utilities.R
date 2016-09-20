library(RPostgreSQL)
drv <- dbDriver("PostgreSQL")
db <- dbConnect(drv, host="mtest.ccsriudvs1y0.us-east-1.rds.amazonaws.com",
                port="5432",
                dbname="metaxcan",
                user="metaxcan_ro",
                password="M3t4xc4n")

db_prs <- dbConnect(drv, host="mtest.ccsriudvs1y0.us-east-1.rds.amazonaws.com",
                    port="5432",
                    dbname="prs",
                    user="metaxcan_ro",
                    password="M3t4xc4n")

db_prs2 <- dbConnect(drv, host="mtest.ccsriudvs1y0.us-east-1.rds.amazonaws.com",
                    port="5432",
                    dbname="prs2",
                    user="metaxcan_ro",
                    password="M3t4xc4n")

db_v6p_hapmap <- dbConnect(drv, host="mtest.ccsriudvs1y0.us-east-1.rds.amazonaws.com",
                           port="5432",
                           dbname="v6p_hapmap",
                           user="metaxcan_ro",
                           password="M3t4xc4n")

## dbListTables(db)
## dbListFields(db,'metaxcan_result')
## dbDisconnect(db)

query_base <- paste0( "SELECT ",
                      " g.gene,",
                      " g.gene_name,",
                      " m.zscore,",
                      " m.effect_size,",
                      " m.pval,",
                      " p.tag as phenotype,",
                      " t.tissue as tissue,",
                      " m.pred_perf_R2,",
                      " m.pred_perf_pval,",
                      " m.n,",
                      " m.model_n",
                      " FROM gene AS g ",
                      " INNER JOIN metaxcan_result AS m ON g.id = m.gene_id ",
                      " INNER JOIN tissue AS t ON t.id = m.tissue_id ",
                      " INNER JOIN pheno AS p ON p.id = m.pheno_id " )

query.pheno.tissue <- function(pheno,tissue)
{
  query <- paste0(query_base,
                  "WHERE p.tag like '", pheno,"%'",
                  " and m.pred_perf_R2 > ", Rthres,
                  "and t.tissue like '" %&% tissue %&% "%'")
  dbGetQuery(db,query)
}

query.pheno <- function(pheno,Rthres=NA,connection=db)
{
  query <- paste0(query_base,
                  "WHERE p.tag like '", pheno, "%'")
  if (!is.na(Rthres)) {
    query <- paste0(query," and m.pred_perf_R2 > ", Rthres)
  }
  dbGetQuery(connection,query)
}

query.all <- function(Rthres=0.)
{
  query <- paste0(query_base,
                  " and m.pred_perf_R2 > ", Rthres)
  dbGetQuery(db,query)
}

build_data <- function(connection=db, phenos = pheno.selected$pheno) {
  d <- data.frame(pred_perf_r2=numeric(), 
                  pred_perf_pval = numeric(),
                  zscore = numeric(),
                  phenotype = character(),
                  tissue = character(),
                  gene = character(),
                  gene_name = character())
  for(phenoname in phenos)
  {
    data <- query.pheno(phenoname, connection=connection)
    append <- data.frame(pred_perf_r2 = data$pred_perf_r2,
                         pred_perf_pval = data$pred_perf_pval,
                         zscore = data$zscore,
                         pval = data$pval,
                         phenotype = data$phenotype,
                         tissue = data$tissue,
                         gene = data$gene,
                         gene_name = data$gene_name,
                         stringsAsFactors=FALSE)
    d <- rbind(d, append)
  }
  return(d)
}

get_all_pheno_names <- function(connection = db) {
  p <- dbGetQuery(connection, "SELECT tag FROM pheno")
  return(p$tag)
}