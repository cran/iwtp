## file: resam.R
## Date created: 2009-11-09
## Date modified: 2010-5-7
## function: iwtp.sample
## input:
## output:
## note: resample, recaluating frequences
iwtp.sample <- function(obj) {

	## define local functions
	
	## function do.sam
	## input: th,  numbers of times of stated intervals
	## output: numbers of times randomly selected for the corresponding intervals
	do.sam <- function(th) {
		tr <- rep.int(1:length(th),times=th)
	
		k <- length(tr)
		# random generating k integers between 1 and k
		rnd <- sort(sample(1:k,k,replace=T))

		ff <- c(tr[rnd])
		##counting
		tab <- table(ff)
		tl <- length(tab)
		count <- rep.int(0,times=length(th))

		ndx <- as.integer(names(tab))
		fr <- as.integer(tab[1:tl])	
		count[ndx] <- fr

		return(count)
	}
	
	th <- obj$intv$n.intv
	new.th <- do.sam(th)
		
  	#new.intv <- obj$intv
  	#new.intv$n.intv <- new.th
  	#print(new.intv)
  
	#out.list <- list(intv=new.intv,freqs.divs=obj$freqs.divs)
	#class(out.list) <-c("iwtp.interv.ties")
	out.list <- list(th=new.th)
	# for debug only
	#out.list <- list(th=th)
	return (out.list)
}
## eof function iwtp.resam


## function: resam
## input:
## output:
## note:
resam <- function(obj,size=10,seed=1234,plot=FALSE,bounds=list(),limits=list()) {
	if(class(obj)!="iwtp")
		stop("obj must be object of class 'iwtp'") 
		
	set.seed(seed)
	n <- size
	run.start <- Sys.time() 
	if (class(obj$sf) == "iwtp.interv.wb") {
		# frequence function parameters
		if(length(limits) < 1L) {
			limits <- defaultPara(mode=obj$mode)		
		}
		
		#x <- lapply(eval(limits), function(p) p[3]) # initial values for parameters
		#names(x) <- names(eval(limits))
		#class(x) <- paste("divs", freq.mode, sep=".")
		#print(x)
		#iwtp.intv <- Dofreq(x, intv, limits)
		
		if(missing(bounds) || length(bounds) < 1L)
			bounds <- list(
					a=c(lower=80,upper=150,init=70),
					b=c(lower=0.5,upper=2.5,init=1.2))
		
		#abs <- matrix(rep(0,2),ncol=2)
		#nparam = ifelse(obj$mode=="ties",2,3)
		nparam <- length(obj$sf$par) + 1 
		abs <- matrix(rep(0,size*nparam),ncol=nparam)
		for(i in 1: n) {
			th.sam <- iwtp.sample(obj$intv)
			ma <- maxer.wb(obj$intv,th=th.sam$th,bounds=bounds)
			abs[i,] <- c(ma$par, ma$value) 
		}

		#abs <- abs[-1,]
		if(length(eval(limits))>0){
			colnames(abs) <- c(c("a","b"),names(eval(limits)),"llike")
		} else {
			colnames(abs) <- c("a","b", "llike")
		}
		## calculate mean wtp 
		m0 <- mean.wb(obj$sf$par)
	
		## meanWTPs for resamplings
		mws <- sort(apply(abs,1,mean.wb))
		devm <- mws - m0 ;

		out <- list(par=abs,sf=obj$sf,size=size,mean.mwtp=mean(mws),
				quantile.mwtp=quantile(mws,probs=c(0.05,0.5,0.95), type=1),dev.mwtp=devm,run.time=Sys.time()-run.start)
		class(out) <- c("iwtp.resam.wb","iwtp.resam")
		
		if(plot)
			plot(out)
		
		return (out)
	}

	if (class(obj$sf) == "iwtp.interv.wemix") {
		# frequence function parameters
		if(length(limits) < 1L) {
			limits <- defaultPara(mode=obj$mode)		
		}

		if(missing(bounds) || length(bounds) < 1L)
			bounds <- list(
					p=c(lower=0.5,upper=0.99,init=0.6),
					a=c(lower=50,upper=300,init=70),
					b=c(lower=1,upper=2.5,init=1.2),
					m=c(lower=200,upper=300,init=210))
		
		#pabms <- matrix(rep(0,4),ncol=4)
		#nparam = 4
		#nparam = ifelse(obj$mode=="ties",4,5)
		nparam <- length(obj$sf$par) + 1 	
		pabms <- matrix(rep(0,size*nparam),ncol=nparam)
		for(i in 1: n) {
			th.sam <- iwtp.sample(obj$intv)
			ma <- maxer.wemix(obj$intv,th=th.sam$th,bounds=bounds)
			pabms[i,] <- c(ma$par, ma$value)
		}

		#pabms <- pabms[-1,]
		#if(obj$mode=="ties")
		if(length(eval(limits))>0) {
			colnames(pabms) <- c(c("p","a","b","m"), names(eval(limits)),"llike")
		} else {
			colnames(pabms) <- c("p","a","b","m","llike")
		}
		## calculate mean wtp
		m0 <- mean.wemix(obj$sf$par)
	
		## meanWTPs
		mws <- sort(apply(pabms,1,mean.wemix))
		devm <- mws - m0 ;

		out <- list(par=pabms,sf=obj$sf,size=size,mean.mwtp=mean(mws),
				quantile.mwtp=quantile(mws,probs=c(0.05, 0.5, 0.95), type=1),dev.mwtp=devm,run.time=Sys.time()-run.start)
		class(out) <- c("iwtp.resam.wemix","iwtp.resam")
		return (out)
	}

}
## eof function resam

plot.iwtp.resam.wb <- function(x,xlab="",ylab="", ...) {
	##if(plot.type =="scatter") {
		abs <- x$par
		## set graph window to square
		par(pty = "s")
		plot(abs,pch=20,col="red",xlab=xlab,ylab=ylab,bty="l",xaxs="i",yaxs="i",...)
	##}
}


plot.iwtp.resam.wemix <- function(x,xlab="",ylab="", ...) {
	##if(plot.type =="scatter") {
		abs <- x$par
		plot(abs[2:3],pch=20,col="red",xlab=xlab,ylab=ylab,xaxs="i",yaxs="i", )
	##}
}


QQplot <- function(x, ...) {UseMethod("QQplot")}
##QQplot.default <- function(x,..) {
##	QQplotStepwise(x,...)
##}


QQplot.iwtp.resam <- function(x,x2=NULL,xlim=NULL,ylim=NULL,xlab="",ylab="",col="red",lty=1,lwd=1,asp=0.8,...) {

	if(missing(x2) || is.null(x2)){
		QQplotStepwise(x$dev.mwtp,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,col=col,lty=lty,lwd=lwd,asp=asp,...)
		}
	else
		QQplot2Stepwise(x$dev.mwtp,x2$dev.mwtp,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,col=col,lty=lty,lwd=lwd,asp=asp,...)	
}


print.iwtp.resam.wb <- function(x, ...) {
	x <- x
	
	cat("Survival function: Weibull\n\n") 

	cat("Resampling size: ")
	cat(x$size) 
	cat("\n")
	
	cat("Mean of mean WTP: ")
	cat(round(x$mean.mwtp,2))

	cat("\n")
	cat("Quantile of mean WTP: \n")
	print(round(x$quantile.mwtp,2))
		
	cat("\n")
	cat("Run time: ")
	cat(format(x$run.time, format = "%X", units = "auto"))
	cat("\n")
}


print.iwtp.resam.wemix <- function(x, ...) {
	
	x <- x
	
	cat("Survival function: mixture of Weibull and exponential\n")
		
	cat("Size of resamplings: ")
	cat(x$size) 
	cat("\n\n")
		
	cat("Mean of mean WTP: ")
	cat(round(x$mean.mwtp,2))

	cat("\n")
	cat("Quantile of mean WTP: \n")
	print(round(x$quantile.mwtp,2))
		
	cat("\n\n")
	cat("Run time: ")
	cat(format(x$run.time, format = "%X", units = "auto"))
	cat("\n\n")
}

