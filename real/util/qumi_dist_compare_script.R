# Compute distances between normalized counts and UMI counts
# for all the different normalizations in the assays of a SingleCellExperiment.
# This script is meant to be run from top level project directory
# eg, Rscript ./real/util/quminorm_script.R ./real/vieira_2019/data/02_sce_qumi.rds ./real/vieira_2019/results/qumi_distance.txt
args<-commandArgs(trailingOnly=TRUE)
ifile<-args[1]
ofile<-args[2]

suppressPackageStartupMessages(library(SingleCellExperiment))
source("./algs/quminorm.R")

sce<-readRDS(ifile)
res<-distance_compare(sce)
write.table(res,file=ofile,quote=FALSE,row.names=FALSE)
