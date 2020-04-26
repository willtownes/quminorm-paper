# functions for fitting and quantile normalizing to poisson-lognormal
#library(slam)
#library(Matrix)
#library(sads)

# rpoilog<-function(n,mu,sig){
#   #same as sads::rpoilog but doesn't involve the quantile function
#   lambda<-rlnorm(n,meanlog=mu,sdlog=sig)
#   rpois(n,lambda)
# }

dpoilog<-function(x,mu,sig,log=FALSE,quadpts=1000){
  #Compute PMF of Poisson-lognormal
  #deprecated, use sads::dpoilog instead!
  #same functionality as poilog::dpoilog
  #mu,sig are same as in poilog::dpoilog
  #if quadpts is a number, the quadrature points are recomputed based on
  #current parameter values
  #if quadpts is a list, we assume it is the result of a call to
  #quadpts<-statmod::gauss.quad.prob(m, dist="normal", mu=0, sigma=1)
  #where m is the number of quadrature points
  #ie that the quadpts are pre-computed based on the Normal(0,1) distribution
  if(is.list(quadpts)){
    q<-quadpts
    lam<-exp(mu+sig*q$nodes)
  } else {
    q<-statmod::gauss.quad.prob(quadpts,dist="normal",mu=mu,sigma=sig)
    lam<-exp(q$nodes)
  }
  w<-log(q$weights)
  f<-function(xi){dpois(xi,lam,log=TRUE)}
  #approximate the log of the integral
  g<-function(xi){matrixStats::logSumExp(w+f(xi))} 
  #log PMF for all data points x and all quadrature points t
  #matrix with nrow=quadpts, ncol=length(x)
  lpmf<-vapply(x,g,FUN.VALUE=0.0)
  if(log==TRUE){ return(lpmf) } else { return(exp(lpmf)) }
}

llcurve_poilog<-function(xmax,lpar,add=TRUE,use_sads=TRUE,quadpts=1000,...){
  #Draw the PMF curve on log-log axes for a Poisson-lognormal distribution 
  #Curve goes from zero to xmax
  #lpar are the mu,sigma parameters (log scale, default for sads and poilog packages)
  if(use_sads){
    f<-function(t){
      sads::dpoilog(floor(expm1(t)),mu=lpar[1],sig=lpar[2],log=TRUE)
    }
  } else {
    f<-function(t){
      dpoilog(floor(expm1(t)),mu=lpar[1],sig=lpar[2],log=TRUE,quadpts=quadpts)
    }
  }
  curve(f,from=0,to=log1p(xmax),add=add,...)
}

llcurve_nb<-function(xmax,lpar,add=TRUE,...){
  #Draw the PMF curve on log-log axes for a negative binomial distribution
  #Curve goes from zero to xmax
  #lpar are the size and mu parameters
  f<-function(t){
    dnbinom(floor(expm1(t)),size=lpar[1],mu=lpar[2],log=TRUE)
  }
  curve(f,from=0,to=log1p(xmax),add=add,...)
}

poilog_mle<-function(x,om="BFGS",...){
  #fit<-poilog::poilogMLE(x,startVals=st,zTrunc=FALSE,method=om,...)
  #mle<-fit$par
  #sads more stable than poilog
  fit<-sads::fitpoilog(x,trunc=NULL,method=om,skip.hessian=TRUE,...)
  mle<-coef(fit) #fit$par
  attr(mle,"loglik")<-as.numeric(logLik(fit)) #fit$logLval
  mle
}

nb_mle<-function(x,...){
  fit<-fitdistrplus::fitdist(x,"nbinom",keepdata=FALSE,...)
  mle<-coef(fit)
  attr(mle,"loglik")<-logLik(fit)
  mle
}

mle_matrix<-function(m,lik=c("poilog","nb"),...){
  #m a matrix with samples in the columns
  #returns a data frame with nrows=ncols(m)
  #result includes MLEs for each model, log likelihood, and BIC
  lik<-match.arg(lik)
  mle_funcs<-list(poilog=poilog_mle,nb=nb_mle)
  f<-mle_funcs[[lik]]
  mle_func<-function(x){
    tryCatch({
      mle<-f(x,...)
      c(mle,loglik=attr(mle,"loglik"))
    },
    error=function(e){
      rep(NA,3)
    })
  }
  if(is(m,"sparseMatrix")){
    apply_func<-function(m){
      m<-slam::as.simple_triplet_matrix(m)
      res<-slam::colapply_simple_triplet_matrix(m,mle_func)
      as.data.frame(do.call(rbind,res))
    }
  } else {
    apply_func<-function(m){ as.data.frame(t(apply(m,2,mle_func))) }
  }
  res<-apply_func(m)
  res$bic<- -2*res$loglik + log(nrow(m))*2
  res
}

poilog_mle_matrix<-function(m,...){
  #mle_matrix for poisson-lognormal
  #result includes mu,sigma params
  mle_matrix(m,"poilog",...)
}

nb_mle_matrix<-function(m,...){
  mle_matrix(m,"nb",...)
}

nb_pzero2mu<-function(lpz,size){
  #Assuming the data follows a negative binomial distribution
  #if we fix the size parameter to a specified value
  #we can infer the mean (mu) parameter from the fraction of zeros
  #lpz is a vector of the log of zero fraction for each cell
  #size can be a scalar or vector
  size*expm1(-lpz/size)
}

poilog_pzero2mu<-function(lpz,sig=2.5,lims=c(-200,200)){
  #Assuming the data follows a Poisson-lognormal
  #if we fix the 'sig' parameter to a specified value
  #we can infer the mu parameter from the fraction of zeros
  #mu != mean of the lognormal, but e^mu is median of lognormal
  #mu can be negative.
  #lpz is a vector of the log of zero fraction for each cell
  #lims are the lower,upper bounds for the mu parameter passed to uniroot
  inner<-function(x,s){
    #x is an element of lpz
    if(is.na(s)){ return(NA) }
    f<-function(mu){
      x-sads::dpoilog(0,mu,sig=s,log=TRUE)
    }
    uniroot(f,lims)$root
  }
  if(length(sig)==1){ #single sig parameter for all cells
    return(vapply(lpz,inner,FUN.VALUE=1.0,s=sig))
  } else { #different tail parameter for each cell
    return(mapply(inner,lpz,sig))
  }
}

################### functions for working with maxima #########################

# note, these functions tend to be numerically unstable for large "n"
# recommend not using them.

poilog_max_cdf<-function(xpts,n,mu,sig,logscale=FALSE){
  #compute the CDF of the max statistic for iid Poisson-lognormals
  #n=number of data points
  #xpts=a grid of values to evaluate the PMF of the max
  #note: n is not the length of xpts!!
  #mu,sig: the parameters of the Poisson-lognormal distribution
  res<- n*sads::ppoilog(xpts,mu=mu,sig=sig,log.p=TRUE)
  if(logscale){ return(res) } else { return(exp(res)) }
}

poilog_max_pmf<-function(xpts,n,mu,sig,logscale=FALSE){
  #compute the PMF of the max statistic
  #n=number of data points
  #xpts=a grid of values to evaluate the PMF of the max
  #note: n is not the length of xpts!!
  #mu,sig: the parameters of the Poisson-lognormal distribution
  l_cdf<-sads::ppoilog(xpts,mu=mu,sig=sig,log.p=TRUE)
  l_pmf<-sads::dpoilog(xpts,mu=mu,sig=sig,log=TRUE)
  res<- log(n)+(n-1)*l_cdf+l_pmf
  if(logscale){ return(res) } else { return(exp(res)) }
}

poilog_max_quantile<-function(q,n,mu,sig,lims=c(0,1000)){
  #compute the theoretical quantile q of the max statistic
  #of a random sample of n data points drawn from
  #a Poisson-lognormal distribution with parameters mu,sig
  #q: the desired quantile to compute (eg 0.5=median)
  xpts<-seq.int(from=lims[1],to=lims[2])
  l_cdf<-poilog_max_cdf(xpts,n,mu,sig,logscale=TRUE)
  i<-max(which(l_cdf<=log(q)))
  if(i==length(xpts)){ stop("reached upper limit without finding quantile") }
  xpts[i]
}

nb_max_cdf<-function(xpts,n,size,mu,logscale=FALSE){
  #compute the CDF of the max statistic for iid negative binomials
  #n=number of data points
  #xpts=a grid of values to evaluate the PMF of the max
  #note: n is not the length of xpts!!
  #size,mu: the parameters of the negative binomial distribution
  res<- n*pnbinom(xpts,size=size,mu=mu,log.p=TRUE)
  if(logscale){ return(res) } else { return(exp(res)) }
}

nb_max_pmf<-function(xpts,n,size,mu,logscale=FALSE){
  #compute the PMF of the max statistic for iid negative binomials
  #n=number of data points
  #xpts=a grid of values to evaluate the PMF of the max
  #note: n is not the length of xpts!!
  #size,mu: the parameters of the negative binomial distribution
  l_cdf<-pnbinom(xpts,size=size,mu=mu,log.p=TRUE)
  l_pmf<-dnbinom(xpts,size=size,mu=mu,log=TRUE)
  res<- log(n)+(n-1)*l_cdf+l_pmf
  if(logscale){ return(res) } else { return(exp(res)) }
}
