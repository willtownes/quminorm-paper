# This script is meant to be run from the top-level directory
# It requires the existence of quminorm_script_config.json in the dataset directory
# specify base path (to the dataset directory) as the first command line argument
# eg, Rscript ./real/util/quminorm_script.R ./real/tung_2016 data/01_sce_all_genes_all_cells.rds data/02_sce_qumi.rds
fp<-file.path
args<-commandArgs(trailingOnly=TRUE)
bp<-args[1]
ifile<-fp(bp,args[2])
ofile<-fp(bp,args[3])

suppressPackageStartupMessages(library(SingleCellExperiment))
source("./util/misc.R") #parcolapply function
source("./algs/quminorm.R")

#shape parameters for plomax, negbinom, and poilog
#tl<-1.0; size<-.1; sig<-2.0; #sig1<-1.0
cfg<-fp(bp,"quminorm_script_config.json")
if(file.exists(cfg)){
  shp<-jsonlite::fromJSON(cfg)
} else {
  stop(paste(cfg,"not found!"))
}
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

for(mod in names(shp)){
  for(s in shp[[mod]]){
    aname<-paste("qumi",mod,s,sep="_")
    if(aname %in% assayNames(sce)){
      message(paste0(aname,": already present in assayNames"))
    } else {
      message(paste0(aname,": computing QUMI counts"))
      assay(sce,aname)<-pfunc(rc,s,lik=mod)
    }
  }
  saveRDS(sce,file=ofile)
}
