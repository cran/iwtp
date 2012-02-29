# file: llikewe.R
# function maxer.wemix
# by wenchao zhou
# Date created: 2010-03-09
# Date modified: 2010-03-19	
# Date modified: 2010-05-07
#####################################################################################
maxer.wemix <- function(x,...) {
	UseMethod("maxer.wemix")
}
## define function maxer.wemix for mixture of Weibull and exponential survial function
maxer.wemix.iwtp.interv <- function(obj,th=obj$intv$n.intv,freq.redo=TRUE,
		bounds=list(
			p=c(lower=0.5,upper=0.99,init=0.6), 
			a=c(lower=70,upper=300,init=75)),
			b=c(lower=1,upper=2.5,init=1.4),
			m=c(lower=200,upper=300,init=210)) {
	## define local functions
	qwemix <- function(v,p,a,b,m) {
		sf.wemix(v[1],p,a,b,m) - sf.wemix(v[2],p,a,b,m)
	}

	llikWE <- function(x, obj=obj,th=obj$intv$n.intv, freq.redo=TRUE) {
		if(class(obj) != "iwtp.interv")
 			stop("Argument `obj' must be class of `iwtp.interv'")

		p <- x[1]
	  	a <- x[2]
  		b <- x[3]
		m <- x[4]

	    ## retrieval data from obj, which is object of class 'iwtp.interv'

		#th <- obj$intv$n.intv
	  	ndx <- obj$intv$ndx.divs
	  	divs <- obj$intv$divs
	  	freqs <- obj$freqs.divs
		ncoef <- length(obj$coef)  # para of freq function
		if(ncoef >0) {
			c <- x[5:(4+ncoef)]
			coef <- list(c)
			names(coef) <- names(eval(obj$freq.limits))
			class(coef) <- paste("divs", obj$freq.mode, sep=".")
	  		freqs <- Dofreq(coef,obj$intv,obj$freq.limits)$freqs.divs
		}
		
		#print(freqs)
	    ## formula for log likelihood
  		sum <- 0
#		for(i in 1:length(th)) {
		for(i in which(th>0)) {
			s <- 0
			for(j in ndx[i,1]:ndx[i,2])  {
			 	s <- s + freqs[i,j]*qwemix(divs[j,],p,a,b,m)
		 	}
		 	#print(s)
		 	if(s >0)	
				sum <- sum + th[i] * log(s)
	    }

  	return (sum)	

	}
	## end local functions

	
	## initialize
	if(length(bounds) <4L)
		stop("The wemix distribtution have four parameters!", call. = FALSE)
	para.names <- c("p", "a", "b", "m")
	names(bounds) <- para.names
	if(length(bounds$p) < 3L || length(bounds$a) < 3L || length(bounds$b) < 3L || length(bounds$m) < 3L) {
		par <- sapply(bounds, function(x) mean(x))  # middle point
	} else
	{
		par <- sapply(bounds,function(x) x[3])
	}
	lower <- sapply(bounds, function(x) x[1])
	upper <- sapply(bounds, function(x) x[2])
	
	
	# append paras from frequence
	if(length(obj$coef)>0){
		coef.lim <- eval(obj$coef.limits)
		par <- c(par,sapply(coef.lim, function(x) x[3]))
		lower <- c(lower, sapply(coef.lim, function(x) x[1]))
		upper <- c(upper, sapply(coef.lim, function(x) x[2]))
	}
	names(par) <- NULL
	names(lower) <- NULL
	names(upper) <- NULL
	# end append
	
## call optimimizer optim 
	max.llikWE <- optim(par,llikWE,obj=obj,th=th,method="L-BFGS-B",
				freq.redo=freq.redo,lower=lower,upper=upper,
				control=list(fnscale=-1,factr=1e-14))	

	out <- list(par=c(max.llikWE$par),value=max.llikWE$value)
	class(out) <- c("iwtp.interv.wemix")
	invisible(out)
}

## survival function in form of mixture of weibull and exponential
sf.wemix <- function(x,p,a,b,m) {
		p*exp(-(x/a)^b) + (1-p)*exp(-x/m)
	}
	
## probability function in form of mixture of weibull and exponential
pf.wemix <- function(x,p,a,b,m) {
	p*(b/a)*((x/a)^(b - 1))*exp[-(x/a)^b] + (1 - p)/m*exp[-x/m]
}

##
qwemix <- function(v,p,a,b,m) {
	sf.wemix(v[1],p,a,b,m) - sf.wemix(v[2],p,a,b,m)
}

## mean for mixture of weibull and exponential
mean.wemix <- function(x) {
	x[1]*x[2]*gamma(1+1/x[3]) + (1-x[1])*x[4]
}

##

rwemix <- function(n=1, par) {
	## par: vector, length(par)=4
	####################
	p <- par[1]
	a <- par[2]
	b <- par[3]
	m <- par[4]
	x <- ifelse(runif(1)<= p, 
			a*exp(log(-log(runif(1)))/b), 
			-m*log(runif(1)))
	return(x)
}


