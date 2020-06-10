# This script is meant to be run from top level project directory
# specify base path (to the dataset directory) as the command line argument and the relative path to the SingleCellExperiment as second argument.
# eg, Rscript ./real/util/mles_qumi.R ./real/klein_2015 data/01_sce_all_genes_all_cells.rds
args<-commandArgs(trailingOnly=TRUE)
bp<-args[1]
sce_file<-args[2]

suppressPackageStartupMessages(library(SingleCellExperiment))
source("./util/misc.R") #parcolapply function
source("./algs/poilog.R")
source("./algs/nblomax.R")
fp<-file.path

sce<-readRDS(fp(bp,sce_file))
m<-counts(sce)
chunksize<-50 #average chunk size=50 cells
odir<-fp(bp,"results")
if(!dir.exists(odir)){dir.create(odir)}

ofile<-fp(odir,"mle_poilog.txt")
if(!file.exists(ofile)){
  message("computing Poisson-lognormal MLEs")
  res<-parcolapply(m,poilog_mle_matrix,chunksize=chunksize)
  write.table(res,file=ofile,quote=FALSE)
}

ofile<-fp(odir,"mle_nb.txt")
if(!file.exists(ofile)){
  message("computing negative binomial MLEs")
  res<-parcolapply(m,nb_mle_matrix,chunksize=chunksize)
  write.table(res,file=ofile,quote=FALSE)
}

ofile<-fp(odir,"mle_poilomax.txt")
if(!file.exists(ofile)){
  message("computing Poisson-Lomax MLEs")
  res<-parcolapply(m,plomax_mle_matrix,chunksize=chunksize,quadpts=1000,maxtry=10)
  write.table(res,file=ofile,quote=FALSE)
}

ofile<-fp(odir,"mle_merged.txt")
if(!file.exists(ofile)){
  message("merging results and discarding cells with NAs")
  nb<-read.table(fp(odir,"mle_nb.txt"),header=TRUE)
  lx<-read.table(fp(odir,"mle_poilomax.txt"),header=TRUE)
  ln<-read.table(fp(odir,"mle_poilog.txt"),header=TRUE)
  good<-!is.na(lx$tail) & !is.na(ln$mu) & !is.na(nb$size)
  good<-intersect(rownames(ln)[good],colnames(sce))
  nb<-nb[good,]; lx<-lx[good,]; ln<-ln[good,]
  ln$method<-"poilog"
  lx$method<-"poilomax"
  nb$method<-"negbinom"
  ln$cell<-rownames(ln); rownames(ln)<-NULL
  lx$cell<-rownames(lx); rownames(lx)<-NULL
  nb$cell<-rownames(nb); rownames(nb)<-NULL
  colnames(ln)[colnames(ln)=="mu"]<-"scale"
  colnames(ln)[colnames(ln)=="sig"]<-"shape"
  colnames(lx)[colnames(lx)=="tail"]<-"shape"
  colnames(nb)[colnames(nb)=="mu"]<-"scale"
  colnames(nb)[colnames(nb)=="size"]<-"shape"
  pd<-rbind(ln,lx,nb)
  write.table(pd,file=ofile,quote=FALSE,row.names=FALSE)
}

ofile<-fp(odir,"compare_maxima.txt")
if(!file.exists(ofile)){
  message("compare simulated maximum from MLEs to true maximum from data")
  mles<-read.table(fp(odir,"mle_merged.txt"),header=TRUE)
  sce<-sce[,as.character(unique(mles$cell))]
  cm<-colData(sce)
  m<-counts(sce)
  true_max<-as.numeric(qlcMatrix::colMax(m))
  rfunc<-function(mthd,scl,shp,n){
    #mthd=method, scl=scale, shp=shape
    #each of the arguments must be individual scalars, not vectors
    #returns the max of a simulated vector of data for the cell
    if(mthd=="poilog"){
      #return(max(sads::rpoilog(n, mu=scl, sig=shp)))
      return(max(poilog::rpoilog(n, mu=scl, sig=shp, keep0=TRUE)))
      #return(max(rpoilog(n,mu=scl,sig=shp)))
    } else if(mthd=="poilomax"){
      return(max(rplomax(n,tail=shp,scale=scl,cut=Inf)))
    } else if(mthd=="negbinom"){
      return(max(rnbinom(n,size=shp,mu=scl)))
    } else {
      stop("invalid method, please choose from poilog, poilomax, or negbinom")
    }
  }
  res<-mapply(rfunc,mles$method,mles$scale,mles$shape,n=nrow(m))
  cpmax<-cbind(mles[,c("method","cell")],sim_max=res,true_max=true_max)
  write.table(cpmax,file=ofile,quote=FALSE,row.names=FALSE)
}
