# This script is meant to be run from top level project directory
# specify base path (to the dataset directory) as the command line argument
# eg, Rscript ./real/util/05_bus2genecount.R ./real/klein_2015
source("./real/util/data_loading.R")
args<-commandArgs(trailingOnly=TRUE)
bus2genecount_all(args[1])
