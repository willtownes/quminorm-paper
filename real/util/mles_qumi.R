# This script is meant to be run from top level project directory
# specify base path (to the dataset directory) as the command line argument
# eg, Rscript ./real/util/mles_qumi.R ./real/klein_2015
args<-commandArgs(trailingOnly=TRUE)
bp<-args[1]

suppressPackageStartupMessages(library(SingleCellExperiment))
source("./util/misc.R") #parcolapply function
source("./algs/poilog.R")
source("./algs/nblomax.R")
fp<-file.path

sce<-readRDS(fp(bp,"data/01_sce_all_genes_all_cells.rds"))
m<-counts(sce)
chunksize<-50 #average chunk size=50 cells

ofile<-fp(bp,"data/mle_poilog.txt")
if(!file.exists(ofile)){
  res<-parcolapply(m,poilog_mle_matrix,chunksize=chunksize)
  write.table(res,file=ofile,quote=FALSE)
}

ofile<-fp(bp,"data/mle_nb.txt")
if(!file.exists(ofile)){
  res<-parcolapply(m,nb_mle_matrix,chunksize=chunksize)
  write.table(res,file=ofile,quote=FALSE)
}

ofile<-fp(bp,"data/mle_poilomax.txt")
if(!file.exists(ofile)){
  res<-parcolapply(m,plomax_mle_matrix,chunksize=chunksize,quadpts=1000,maxtry=10)
  write.table(res,file=ofile,quote=FALSE)
}
