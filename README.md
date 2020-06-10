# Quantile normalization of single-cell RNA-seq read counts without unique molecular identifiers

[![DOI](https://zenodo.org/badge/217110308.svg)](https://zenodo.org/badge/latestdoi/217110308)

This repository contains supporting code to facilitate reproducible analysis. For details see the [biorxiv preprint](https://www.biorxiv.org/content/10.1101/817031v1). If you find bugs please create a github issue. 

The quasi-UMI normalization method is under development as a [standalone R package](https://github.com/willtownes/quminorm).

The purpose of this method is to remove PCR distortion from scRNA-seq read counts by normalizing to quasi-UMIs (QUMIs). QUMIs approximate the true (unmeasured) UMI counts. Once read counts are transformed to QUMIs, the count matrix can be passed to UMI-specific methods for [feature selection and dimension reduction](https://github.com/willtownes/scrna2019).

### Authors

Will Townes and Rafael Irizarry

## Description of Repository Contents

### algs

Implementation of quasi-UMI normalization methods, auxiliary functions, and methods used in comparisons.
* *existing.R* - wrapper functions for PCA, tSNE, ZINB-WAVE, etc.
* *nblomax.R* - density, visualization, and maximum likelihood estimation of compound negative binomial- Lomax families. These are discrete distributions with heavy power-law tails.
* *poilog.R* - density, visualization, and maximum likelihood estimation of compound Poisson-lognormal families. These are discrete distributions with heavy tails (but not as heavy as power law).
* *quminorm.R* - quasi-UMI normalization. The most important function is *quminorm_matrix* which applies the normalization to a matrix of read counts with features in rows, samples in columns. Also of interest is *lpmf* which visualizes heavy-tailed data with log-log plots (histograms with logarithmic binning).

### qumi

Code for making graphs that involved multiple datasets.

### real

Analysis of various real scRNA-seq datasets. The Rmarkdown files can be used to produce figures in the manuscript. The *util* subfolder contains numerous scripts for generating count matrices from FASTQ files, working with [BUS files](https://bustools.github.io/), and running QUMI normalization on large matrices using parallelization on a high performance computing cluster.

### util

Miscellaneous utility functions. 
