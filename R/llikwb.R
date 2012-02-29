# file: llikwb.R
# maxer.wb
# by wenchao zhou
# created: 2009-03-09
# modified: 2010-05-07
#
###################################################################################
maxer.wb <- function(x,...) {
	UseMethod("maxer.wb")
}
## function maxer for Weibull survival function
maxer.wb.iwtp.interv <- function(obj,th=obj$intv$n.intv,freq.redo=TRUE,
		bounds=list(
			a=c(lower=70,upper=150,init=75),
			b=c(lower=1,upper=2.5,init=1.4))){
	##define local functions
	qw <- function(a,b,hl,hr) { 
		sf.wb(hl,a,b) - sf.wb(hr,a,b)
	}

	llikWB <- function(x,obj=obj,th=obj$intv$n.intv,freq.redo=TRUE) {
		if(class(obj) != "iwtp.interv")
 			stop("Argument `obj' must be object of class `iwtp.interv'")

		#print(obj)
	  	a <- x[1]
  		b <- x[2]
  	
	  	## retrieval freqs from obj, which is object of class 'iwtp.interv'
	  	#th <- obj$intv$n.intv
		ndx <- obj$intv$ndx.divs
	  	divs <- obj$intv$divs
	  	freqs <- obj$freqs.divs
		ncoef <- length(obj$coef)
	  	if(ncoef>0) {
	  		c <- x[3:(2+ncoef)]
	  		coef <- list(c)
			names(coef) <- names(eval(obj$freq.limits))
			class(coef) <- paste("divs", obj$freq.mode, sep=".")
	  		freqs <- Dofreq(coef,obj$intv,obj$freq.limits)$freqs.divs
	  	}	
	  
	  	#print(th)
	  	#print(ndx)
	  	#print(divs)
	  	#print(freqs)
	  
	  	## formula for log likelihood
  		sum <- 0
#		for(i in 1:length(th)) {
		for(i in which(th>0)) {
			s <- 0
			for(j in ndx[i,1]:ndx[i,2])  {
			 	  	s <- s + freqs[i,j]*(qw(a,b,divs[j,1],divs[j,2]))
		 		}
		 	if(s >0)	
				sum <- sum + th[i] * log(s)
	    }

		return (sum)
	}

## eof local function

## call optimimizer optim 
    if(length(bounds) <2L)
		stop("The weibull distribtution have two parameters!", call. = FALSE)
	para.names <- c("a", "b")
	names(bounds) <- para.names
	if(length(bounds$a) < 3L || length(bounds$b < 3L)) {
	   	par <- c(mean(bounds$a), mean(bounds$b))
	} else
	{
		  par <-sapply(bounds, function(x) x[3])
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
		
    #browser()
	max.llikWB <- optim(par,llikWB,obj=obj,th=th,method="L-BFGS-B",
			freq.redo=freq.redo,lower=lower,upper=upper,
			control=list(fnscale=-1,factr=1e-14))
	out <- list(par=c(max.llikWB$par),value=max.llikWB$value)
	class(out) <- c("iwtp.interv.wb")
	invisible(out)
}
## eof fucntion maxer.wb

 
#
## survival function in form of weibull
sf.wb <- function(x,a,b) {
	exp(-(x/a)^b)
}

#
## mean for weibull
mean.wb <- function(x) {
	return (x[1]*gamma(1+1/x[2]))
}
