# This script is meant to be run from top level project directory
# specify base path (to the dataset directory) as the command line argument
# eg, Rscript ./real/util/quminorm_script.R ./real/tung_2016/data/01_sce_all_genes_all_cells.rds ./real/tung_2016/data/02_sce_qumi.rds
args<-commandArgs(trailingOnly=TRUE)
ifile<-args[1]
ofile<-args[2]

suppressPackageStartupMessages(library(SingleCellExperiment))
source("./util/misc.R") #parcolapply function
source("./algs/quminorm.R")
fp<-file.path

#shape parameters for plomax, negbinom, and poilog
tl<-1.0; size<-.1; sig<-2.0; #sig1<-1.0
chunksize<-100
cores<-max(parallel::detectCores()-1,1)

sce0<-readRDS(ifile)
if(file.exists(ofile)){
  sce<-readRDS(ofile)
  stopifnot(all(assayNames(sce0) %in% assayNames(sce)))
} else {
  sce<-sce0
}
rc<-assay(sce,"read_counts")
if(!("cpm" %in% assayNames(sce))){
  assay(sce,"cpm")<-rc2cpm(rc)
}
if(!("scran_counts" %in% assayNames(sce))){
  assay(sce,"scran_counts")<-scran_normalize(rc)
}
if(!("census_counts" %in% assayNames(sce))){
  message("computing census counts")
  sce<-census_normalize(sce,assay_name="cpm",cores=cores)
  saveRDS(sce,file=ofile)
}
pfunc<-function(rc,shape,lik,...){
  parcolapply(rc,quminorm_matrix,transpose=FALSE,chunksize=chunksize,cores=cores,shape=shape,lik=lik)
}
if(!("qumi_nb" %in% assayNames(sce))){
  message("computing negative binomial QUMIs")
  #assay(sce,"qumi_nb")<-quminorm_matrix(rc,size,lik="nb")
  assay(sce,"qumi_nb")<-pfunc(rc,size,lik="nb")
}
if(!("qumi_poilog" %in% assayNames(sce))){
  message("computing Poisson-lognormal QUMIs")
  #assay(sce,"qumi_poilog")<-quminorm_matrix(rc,sig,lik="poilog")
  assay(sce,"qumi_poilog")<-pfunc(rc,sig,lik="poilog")
  saveRDS(sce,file=ofile)
}
# if(!("qumi_poilog1" %in% assayNames(sce))){
#   message("computing Poisson-lognormal QUMIs, secondary sigma value")
#   assay(sce,"qumi_poilog1")<-pfunc(rc,sig1,lik="poilog")
#   saveRDS(sce,file=ofile)
# }
if(!("qumi_plomax" %in% assayNames(sce))){
  message("computing Poisson-Lomax QUMIs")
  #assay(sce,"qumi_plomax")<-quminorm_matrix(rc,tl,lik="plomax")
  assay(sce,"qumi_plomax")<-pfunc(rc,tl,lik="plomax")
}
saveRDS(sce,file=ofile)
