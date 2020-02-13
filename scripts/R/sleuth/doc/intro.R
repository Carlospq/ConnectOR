## ----eval=FALSE----------------------------------------------------------
## source("http://bioconductor.org/biocLite.R")
## biocLite("rhdf5")

## ----eval=FALSE----------------------------------------------------------
## install.packages("devtools")

## ----eval=FALSE----------------------------------------------------------
## devtools::install_github("pachterlab/sleuth")

## ------------------------------------------------------------------------
library("sleuth")

## ------------------------------------------------------------------------
base_dir <- "~/Downloads/cuffdiff2_data_kallisto_results"

## ------------------------------------------------------------------------
sample_id <- dir(file.path(base_dir,"results"))

## ------------------------------------------------------------------------
sample_id

## ------------------------------------------------------------------------
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results", id, "kallisto"))
kal_dirs

## ------------------------------------------------------------------------
s2c <- read.table(file.path(base_dir, "hiseq_info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c

## ------------------------------------------------------------------------
s2c <- dplyr::mutate(s2c, path = kal_dirs)

## ------------------------------------------------------------------------
print(s2c)

## ----eval=TRUE-----------------------------------------------------------
so <- sleuth_prep(s2c, ~ condition)

## ----eval=TRUE-----------------------------------------------------------
so <- sleuth_fit(so)

## ------------------------------------------------------------------------
so <- sleuth_fit(so, ~1, 'reduced')

## ----eval=TRUE-----------------------------------------------------------
so <- sleuth_lrt(so, 'reduced', 'full')

## ----eval=TRUE-----------------------------------------------------------
models(so)

## ----eval=FALSE----------------------------------------------------------
## source("http://bioconductor.org/biocLite.R")
## biocLite("biomaRt")

## ----eval=TRUE-----------------------------------------------------------
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'ensembl.org')

## ---- eval=TRUE----------------------------------------------------------
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

## ---- eval=FALSE---------------------------------------------------------
## sleuth_live(so)

## ------------------------------------------------------------------------
results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')

## ------------------------------------------------------------------------
so <- sleuth_prep(s2c, ~condition, target_mapping = t2g,
  aggregation_column = 'ens_gene')

## ----eval=FALSE----------------------------------------------------------
## # set the number of available cores to 4
## options(mc.cores = 4L)

