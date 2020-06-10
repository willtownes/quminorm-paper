#Quasi-UMIs:
#estimate scale parameter using only the fraction of zeros
#quantile normalization to a Poisson-Lomax or Geometric-Lomax
source("./algs/nblomax.R")
source("./algs/poilog.R")

make_cdf_nz<-function(thresh,dfunc,maxval=1e6){
  #dfunc=some log-pmf function accepting a single argument (the data)
  #let cdf be the cdf corresponding to dfunc
  #thresh is a probability value, where if cdf>thresh,...
  #...it implies 1-cdf_nz is less than 1/number of nonzero data points
  #beyond this point the quantile function resolves to less than one data point
  #so no point in computing cdf above threshold
  #VALUE: cdf_nz: the zero-truncated version of cdf
  #cdf_nz can then be used to quantile normalize data to the distribution of dfunc
  lthresh<-log(thresh)
  lo<-0; hi<-100
  lcdf<-NULL
  lcdf_max<- -Inf
  ctr<-0 #counter
  while(lcdf_max<=lthresh && hi<maxval){
    ctr<-ctr+1
    z<-seq.int(lo,hi-1) #0:99, 100:199, 200:399, ...
    lpmf<-dfunc(z)
    #update the first pmf value by adding the max cdf to it
    lpmf[1]<-logspace_add(lpmf[1],lcdf_max)
    lcdf<-c(lcdf,log_cumsum_exp(lpmf))
    lcdf_tail<-tail(lcdf,2) #make sure not to change the "2", it is crucial!
    if(diff(lcdf_tail)==0){
      stop("CDF is not increasing, PMF function may have numerical problems!")
    }
    lcdf_max<-lcdf_tail[2] #lcdf_tail must be length 2!!
    lo<-hi #100
    hi<-2*hi #200
  }
  if(hi>maxval){
    stop("Exceeded max value, PMF function may have numerical problems!")
  }
  #remove extra elements of cdf that extend beyond the threshold
  k<-min(which(lcdf>lthresh,arr.ind=TRUE))
  #renormalize CDF for nonzero values only
  #cdf[1] is P(x=0)
  pmf0<-exp(lcdf[1])
  (exp(lcdf[2:k])-pmf0)/(1-pmf0)
}

quminorm_inner<-function(xnz,cdf_nz,nnz=length(xnz)){
  #xnz a vector of positive integers
  #cdf_nz a CDF for some zero-truncated distribution
  #value: the quantile normalized version of xnz
  rk_th<-ceiling(cdf_nz*nnz) #theoretical ranks, increasing order
  dups<-duplicated(rk_th) #TRUE occurs for first duplicate, then all FALSE
  targets<-seq_along(rk_th)[!dups] #the actual qumi values we map to.
  rk_th<-rk_th[!dups]
  xrk<-rank(xnz,ties.method="min") #convert data to empirical ranks
  #xmap gives indices within targets where the quantile function points to
  xmap<-findInterval(xrk,rk_th,left.open=TRUE)+1 #intervals are (lo,hi]
  targets[xmap]
}

quminorm_poilog<-function(x,shape,sc=NULL,quadpts=1000,err2na=TRUE){
  #shape=sig, sc=mu
  #sig,mu are params of lognormal as in poilog::dpoilog and sads::dpoilog
  #returns a quantile normalized version of the data x with zeros unchanged
  n<-length(x)
  nzi<-which(x>0)
  xnz<-x[nzi]
  nnz<-length(xnz)
  if(is.null(sc)){
    #assumes at least one gene is a zero count
    lpz<-log(n-nnz)-log(n) #log(fraction of zeros)
    sc<-poilog_pzero2mu(lpz,sig=shape)
  }
  dfunc<-function(x){ dpoilog(x,mu=sc,sig=shape,quadpts=quadpts,log=TRUE) }
  pmf0<-exp(dfunc(0))
  #threshold for cdf on the regular scale such that 
  #zero truncated cdf (cdf_nz) extends beyond the 1-1/nnz threshold
  thresh<-(1-1/nnz)*(1-pmf0)+pmf0
  if(err2na){
    cdf_nz<-tryCatch(make_cdf_nz(thresh,dfunc),error=function(e){NULL})
  } else {
    cdf_nz<-make_cdf_nz(thresh,dfunc)
  }
  if(is.null(cdf_nz)){
    x[nzi]<-NA
  } else {
    x[nzi]<-quminorm_inner(xnz,cdf_nz,nnz)
  }
  x #quantile normalized version of x
}

quminorm_plomax<-function(x,shape,sc=NULL,quadpts=1000){
  #shape=tail, sc=scale
  #tail,sc are the lomax (power law) tail and scale params, see nblomax.R
  #returns a quantile normalized version of the data x with zeros unchanged
  n<-length(x) 
  nzi<-which(x>0)
  xnz<-x[nzi]
  nnz<-length(xnz)
  if(is.null(sc)){
    #assumes at least one gene is a zero count
    lpz<-log(n-nnz)-log(n) #log(fraction of zeros)
    sc<-pglomax_pzero2scale(lpz,lik="poisson",tail=shape,quadpts=quadpts)
  }
  dfunc<-function(x){ dplomax(x,tail=shape,scale=sc,quadpts=quadpts,log=TRUE) }
  pmf0<-exp(dfunc(0)) 
  #threshold for cdf on the regular scale such that 
  #zero truncated cdf (cdf_nz) extends beyond the 1-1/nnz threshold
  thresh<-(1-1/nnz)*(1-pmf0)+pmf0
  cdf_nz<-make_cdf_nz(thresh,dfunc)
  x[nzi]<-quminorm_inner(xnz,cdf_nz,nnz)
  x #quantile normalized version of x
}

quminorm_nb<-function(x,shape,sc=NULL,quadpts=NULL){
  #shape=size, sc=mu
  #size,mu are params of negative binomial as in dnbinom
  #returns a quantile normalized version of the data x with zeros unchanged
  #quadpts is ignored, included only for consistency with other quminorm funcs
  n<-length(x)
  nzi<-which(x>0)
  xnz<-x[nzi]
  nnz<-length(xnz)
  if(is.null(sc)){
    #assumes at least one gene is a zero count
    lpz<-log(n-nnz)-log(n) #log(fraction of zeros)
    sc<-nb_pzero2mu(lpz,size=shape)
  }
  #dfunc<-function(x){ dnbinom(x,size=shape,mu=sc,log=TRUE) }
  #pmf0<-exp(dfunc(0)) 
  pmf0<-dnbinom(0,size=shape,mu=sc,log=FALSE)
  #threshold for cdf on the regular scale such that 
  #zero truncated cdf (cdf_nz) extends beyond the 1-1/nnz threshold
  thresh<-(1-1/nnz)*(1-pmf0)+pmf0
  qstop<-qnbinom(thresh,size=shape,mu=sc) #an integer
  #note after normalization all values should be strictly less than qstop
  cdf<-pnbinom(seq_len(qstop),size=shape,mu=sc)
  #renormalize CDF for nonzero values only
  cdf_nz<-(cdf-pmf0)/(1-pmf0)
  x[nzi]<-quminorm_inner(xnz,cdf_nz,nnz)
  x #quantile normalized version of x
}

quminorm_matrix<-function(m,shape,lik=c("poilog","plomax","nb"),quadpts=1000){
  #m a matrix with samples in the columns
  #returns a matrix of same dims that has been quantile normalized
  #shape should be a scalar
  #... additional args such as quadpts passed to quminorm functions
  lik<-match.arg(lik)
  qfunc<-switch(lik,poilog=quminorm_poilog,plomax=quminorm_plomax,nb=quminorm_nb)
  res<-0*m
  for(i in seq_len(ncol(m))){
    res[,i]<-qfunc(m[,i],shape,quadpts=quadpts)
  }
  res
}

scran_normalize<-function(rc){
  #rc a read counts matrix or sparse Matrix
  cl<-scran::quickCluster(rc)
  sz<-scran::calculateSumFactors(rc,clusters=cl)
  t(t(rc)/sz)
}

rc2cpm<-function(rc){
  t(t(rc)/Matrix::colSums(rc))*1e6
}

census_normalize<-function(sce,assay_name="cpm",...){
  #sce a single cell experiment object
  #returns the same sce object but with new assay "census_counts"
  # First create a CellDataSet from the relative expression levels
  cds<-scran::convertTo(sce,"monocle",assay.type=assay_name)
  assay(sce,"census_counts")<-monocle::relative2abs(cds,method="num_genes",...)
  sce
}

#this commented out function is deprecated in favor of the one below
#I tested to make sure they give the same result
#the function below this one is much faster!
# distance_compare<-function(sce0,baseline="counts"){
#   #compare the distances between all assays in sce0
#   #versus the baseline assay
#   #returns a data frame with nrows=ncol(sce0)
#   f<-function(j,baseline="counts"){
#     sce<-sce0[,j]
#     nz<-which(drop(assay(sce,baseline)>0))
#     sce<-sce[nz,]
#     g<-function(a){log(drop(assay(sce,a)))}
#     dat<-vapply(assayNames(sce),g,FUN.VALUE=rep(0.0,length(nz)))
#     as.matrix(dist(t(dat)))[,baseline]
#   }
#   res<-vapply(seq_len(ncol(sce0)),f,FUN.VALUE=rep(0.0,length(assayNames(sce0))))
#   res<-as.data.frame(t(res))
#   res$cell<-colnames(sce0)
#   res
# }

distance_compare<-function(sce,baseline="counts"){
  umi<-assay(sce,baseline)
  sparse0<-is(umi,"dgCMatrix")
  if(sparse0){
    umi<-Matrix::drop0(umi)
    z<- umi>0
    umi@x<-log(umi@x)
  } else {
    z<- umi>0
    umi<-log(umi)
  }
  an<-setdiff(assayNames(sce),baseline)
  res<-matrix(0,nrow=ncol(sce),ncol=length(an))
  colnames(res)<-an
  res<-cbind(as.data.frame(res),cell=colnames(sce))
  for(a in an){
    m<-z*assay(sce,a) #make sure sparsity patterns agree
    if(sparse0 && is(m,"dgCMatrix")){
      m<-Matrix::drop0(m)
      m@x<-log(m@x)
      res[,a]<-sqrt(Matrix::colSums((m-umi)^2))
    } else {
      m<-log(m)
      res[,a]<-sqrt(colSums((m-umi)^2,na.rm=TRUE))
    }
  }
  res
}

# format_qumi_assayname<-function(mod,shp){
#   #given a QUMI model and a shape parameter value,
#   #return a string that can be used as an assay name for storing the
#   #QUMI normalized counts in a SingleCellExperiment object.
#   #mod=a string like "nb", "poilog", or "plomax"
#   #shp= the value of the shape parameter, a number
#   #returns a string like nb-0_1 (if shape=0.1)
#   #or poilog-2_4 (if shape=2.4)
#   #or plomax-1 (if shape=1.0)
#   nm<-paste0(mod,shp,sep="-")
#   sub(".","_",nm,fixed=TRUE)
# }

################ visualization functions ################

lpmf<-function(x,bw=.25,logbase=NULL,discrete=TRUE,midpt=FALSE,bin_thresh=10,add=FALSE,doplot=TRUE,...){
  #x a vector of counts or continuous values
  #bw the bandwidth of the histogram
  #if midpt=TRUE uses midpoints of bins by geometric mean
  #if midpt=FALSE uses the left side of the bin
  #probs: if discrete, convert x to probabilities by dividing by total, handles exact zeros
  if(is.null(logbase)){ #defaults to base e
    la<-1
    xlab<-"log(1+x)"
  } else {
    la<-log(logbase)
    xlab<-paste0("log",logbase,"(1+x)")
  }
  bw<-bw*la
  if(discrete){
    xt<-table(x)
    u<-as.integer(names(xt))
    xd<-diff(xt)
    #find point where need to start binning
    if(xd[1]>=0){
      #account for possible initial increase in the PMF
      j<-which.max(xt)
      #find point after the maximum at which it drops below bin_thresh counts
      i<-min(which(xt[j:length(xt)]>=bin_thresh))+(j-1)
    } else {
      #switch point from decreasing to increasing
      i<-min(which(xd>=0)) 
    }
  } else {
    i<-1
  }
  if(!is.infinite(i)){
    bmin<-if(i>1){ log(u[i]) } else { log(min(x)) }
    bmax<-log1p(max(x))
    blen<-ceiling((bmax-bmin)/bw)
    bseq<-seq(from=bmin,to=bmax,length.out=max(2,blen))
    breaks<-exp(bseq)
    if(i>1) breaks<-c(u[1:(i-1)],breaks)
  } else {
    breaks<-c(u[-length(u)],exp(log1p(max(u))))
  }
  h<-hist(x,breaks=breaks,right=FALSE,plot=FALSE)
  if(midpt){
    #use geometric mean to get midpoints
    gmean<-function(t){floor(sqrt(breaks[t]*breaks[t+1]))}
    xpts<-vapply(seq_along(breaks[-length(breaks)]),gmean,FUN.VALUE=1.0)
  } else {
    xpts<-breaks[-length(breaks)]
  }
  good<-h$density>0
  xvals<-log1p(xpts[good])/la
  yvals<-log(h$density[good])
  if(doplot){
    if(add){
      lines(xvals,yvals,...)
    } else {
      plot(xvals,yvals,xlab=xlab,ylab="log(density)",...)
    }
  } else {
    return(data.frame(x=xvals,y=yvals))
  }
}

lpmf_xtra<-function(lpmf_obj,connect=TRUE,logbase=NULL,...){
  #lpmf_obj is data frame result of lpmf function with doplot=FALSE
  #makes a fancier version of lpmf with no pseudocount, logarithmic x-axis
  #... args passed to plot function
  res<-lpmf_obj
  if(is.null(logbase)){
    res$x<-expm1(res$x)
  } else {
    res$x<-logbase^res$x - 1
  }
  res$lx<-log(res$x)
  res$lx[1]<- -log(2)
  plot(res$lx,res$y,xaxt="n",...)
  if(connect){
    lines(res$lx[-1],res$y[-1])
  }
  axis(1, at=res$lx, labels=round(res$x,1))
  plotrix::axis.break(1,-log(2)/2)
  plotrix::axis.break(3,-log(2)/2)
}

lsurv<-function(x,...){
  #x is a random sample from some discrete distribution
  #plots the log of empirical survival function vs log1p(unique values)
  xt<-table(x)
  G<-length(x)
  ly<-log(G-cumsum(xt))-log(G)
  ylab<-"log(1-CDF)"
  lx<-log1p(as.integer(names(ly)))
  plot(lx,ly,xlab="log(1+x)",ylab=ylab,...)
}