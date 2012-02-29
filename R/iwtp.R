## file: iwtp.R
## function: iwtp
## Date created: 2009-2-15
## Date modified: 2010-3-23
## Date modified: 2010-05-20
## input: data -- cases of self-selected intervals
#  		  dist --  survival function:  weibull or wemix
##        bm.type -- behaviour model
## output: object of class 'iwtp' 
##
iwtp <- function(data,dist="weibull",bm.type=5,plot=FALSE,bounds=list(),limits=list()) {
	# data check		
	if (is.null(data)) 
	  	stop("Data must be specified!", call. = FALSE)

	if (!is.data.frame(data)) {
    		stop("Data must be a dataframe!", call. = FALSE)
	}
	
	if(dim(data)[2]<3)
		stop("The input data should have at least 3 columns!",call.=FALSE)
	
	if (sum(is.na(data)) > 0) 
      	stop("Missing values exist in the data!", call.=FALSE)

	# check parameter: survival function
	if(!(tolower(dist) %in% c("weibull","wemix")))
		stop(sprintf("survival function '%s' not yet implemented", dist), call. = FALSE)

	# check parameter: behaviour model		
	bms <- paste("BM", 1:5, sep="") # 5 behaviour models implemented
	if(is.numeric(bm.type) && bm.type %in% c(1:5))
		bm.type <- bms[bm.type]
		
	if(!(toupper(bm.type) %in% bms))
		stop(sprintf("modeling behavior type '%s' not yet implemented", bm.type), call.=FALSE)

	freq.mode = as.character(bm.type)

	# set behaviour model function parameters
	if(length(limits) < 1L) {
	   limits <- defaultPara(mode=freq.mode)		
	}
	
	x <- lapply(eval(limits), function(p) p[3]) # initial values for parameters
	names(x) <- names(eval(limits))
	class(x) <- paste("divs", freq.mode, sep=".")
	#print(x)

	data <- as.matrix(data)
	tab <- table(data[,2],data[,3])
	data.mat <- tab2matrix(tab)
	data.df <- as.data.frame(data.mat)
	names(data.df) <- c("interv no","lower value", "upper value","th")
	##
	intv <- interv(data.df)
	iwtp.intv <- Dofreq(x, intv, limits)
    #
	th <- iwtp.intv$intv$n.intv  # number of a given interval being selected
	if(dist == "weibull") {
		if(missing(bounds) || length(bounds) < 1L)
		{
			
			#bounds <- list(
			#		a=c(lower=80,upper=150,init=70),
			#		b=c(lower=0.5,upper=2.5,init=1.2))
			a.lower <- min(data.df[,'upper value'])
			a.upper <- max(data.df[,'upper value'])
			bounds <- list(
					a=c(lower=a.lower,upper=a.upper,init=mean(data.df[,'upper value'])),
					b=c(lower=0.5,upper=20,init=1.2))
			
		}
		sf.par <- maxer.wb(iwtp.intv,th=th,freq.redo=T,bounds=bounds)
	}
	else {
		if(missing(bounds) || length(bounds) < 1L)
		{
			#bounds <- list(
			#		p=c(lower=0.5,upper=0.99,init=0.6),
			#		a=c(lower=50,upper=300,init=70),
			#		b=c(lower=1,upper=2.5,init=1.2),
			#		m=c(lower=200,upper=300,init=210))
			a.lower <- min(data.df[,'upper value'])
			a.upper <- max(data.df[,'upper value'])
			bounds <- list(
				p=c(lower=0.5,upper=0.99,init=0.6),
				a=c(lower=a.lower,upper=a.upper,init=mean(data.df[,'upper value'])),
				b=c(lower=1,upper=10,init=1.2),
				m=c(lower=200,upper=300,init=210))
		}
		sf.par <- maxer.wemix(iwtp.intv,th=th,freq.redo=T,bounds=bounds)
	}
	## empirical survival function, nonparametric estimate
	esf.nonpar <- iwtp.esf(data.df[,3:4])	
	
	# recalculte frequence using optimal coef 
	if(length(x) >0 ) {
		ncoef <- length(x)
		n <- length(sf.par$par)
		coef.optim <- sf.par$par[(n-ncoef+1):n]
		#print(coef.optim)
		x1 <- lapply(coef.optim, function(y) y)
		names(x1) <- names(x)
		class(x1)=paste("divs", freq.mode, sep=".")
		iwtp.intv <- Dofreq(x1,intv,limits)
	}

	# the left value cann't be 0
	mean.turnbull.weibull <- turnbull.weibull(data)
	#
	surv.turnbull <- turnbull(cbind(data[,2],data[,3]))
	mean.turnbull <- auc.icsurv(surv.turnbull)
	#print(mean.me)
	# output
	out <- list(data=data.df,sf=sf.par,esf=esf.nonpar,intv=iwtp.intv,mode=freq.mode,mean.turnbull.weibull=mean.turnbull.weibull,mean.turnbull=mean.turnbull)
	class(out) <- "iwtp"
	
	if(plot)
		esfPlot(out,xlab="",ylab="")
	
	return(out)
}

## eof function iwtp

print.iwtp <- function (x, ...) {
	
	x <- x

	print(x$intv)
	
	if(class(x$sf) == "iwtp.interv.wb") {
		cat("Survival function: Weibull\n") 
		cat("Behavior model:", x$mode, "\n")
		cat("Parameter estimates:\n")
		#print(x$sf$par)
		cat(sprintf("    Scale (a): %f\n", x$sf$par[1]))
		cat(sprintf("    Shape (b): %f\n", x$sf$par[2]))
		if(length(x$intv$coef) >0 ) {
			for(i in 1: length(x$intv$coef))
				cat(sprintf("    Coefficient %s: %f\n", names(x$intv$coef)[i], x$sf$par[i+2]))
		}
		cat(sprintf("Maximum llik: %f\n", x$sf$value))
		cat(sprintf("Mean WTP: %f\n", mean.wb(x$sf$par)))
		cat("\n")
	}


	if(class(x$sf) == "iwtp.interv.wemix") {
		cat("Survival function: the mixture of Weibull/exponential\n") 
		cat("Bahavior model:", x$mode, "\n")
		cat("Parameter estimates:\n")
		#print(x$sf$par)
		cat(sprintf("    mixing parameter p: %f\n", x$sf$par[1]))
		cat(sprintf("    scale parameter a: %f\n", x$sf$par[2]))
		cat(sprintf("    shape parameter b: %f\n", x$sf$par[3]))
		cat(sprintf("    parameter m: %f\n", x$sf$par[4]))
		if(length(x$intv$coef) >0) {
			for(i in 1: length(x$intv$coef))
				cat(sprintf("    Coefficient %s: %f\n", names(x$intv$coef)[i], x$sf$par[i+4]))			
		}
		cat(sprintf("Maximum llik: %f\n", x$sf$value))
		cat(sprintf("Mean WTP: %f\n",mean.wemix(x$sf$par)))
		cat("\n")
	}

	#cat("\n")
	cat("Nonparametric estimate of the mean value of the empirical survival function\n")
	cat("(based on the right-end values of self-selected intervals): ")
	cat(x$esf$mean)
	cat("\n")
	#
	cat("\n")
	cat("Nonparametric estimate of the mean value based on Turnbull estimator\n")
	cat("[see Turnbull, 1976]: ")
	cat(x$'mean.turnbull')
	cat("\n")
	#
	cat("\n")
	cat("Estimated Mean WTP based on parametric Weibull model (taken self-selected\n")
	cat("intervals as Turnbull intervals [see Turnbull, 1976]): ")
	cat(x$'mean.turnbull.weibull')
	cat("\n")

}


## define generic function: esfPlot
esfPlot <- function(x,...) {UseMethod("esfPlot")}

esfPlot.iwtp <- function(x,x2=NULL,xlim=NULL,xlab="",col="red",lty=1,lwd=1,asp=0.6,...) {
	if(missing(x2) || is.null(x2)) {
		iwtp.esfPlot(x$esf$esf,xlim=xlim,xlab=xlab,col=col,lty=lty,lwd=lwd,asp=asp,...)
	}
	else
	{
		iwtp.esf2Plot(x$esf$esf,x2$esf$esf,xlim=xlim,xlab=xlab,col=col,lty=lty,lwd=lwd,asp=asp,...)
	}
}


esfPlot.iwtp.esf <- function(x,x2=NULL,xlim=NULL,xlab="",col="red",lty=1,lwd=1,asp=0.6,...) {
	if(missing(x2) || is.null(x2)) {
		iwtp.esfPlot(x$esf,xlim=xlim,xlab=xlab,col=col,lty=lty,lwd=lwd,asp=asp,...)
	}
	else
	{
		iwtp.esf2Plot(x$esf,x2$esf,xlim=xlim,xlab=xlab,col=col,lty=lty,lwd=lwd,asp=asp,...)
	}
}

iwtp.esfPlot <- function(esf,xlim=NULL,xlab="",col="red",lty=1,lwd=1,asp=0.6,...) {
	if(missing(xlim)||is.null(xlim))
		xlim <- c(0,max(esf[,1]))

	if(missing(xlab)||is.null(xlab))
		xlab = ""
	#if(missing(ylab)||is.null(ylab))
	#	ylab = ""
	## set graph window to square
	par(pty = "s")
	plot.new()
	plot.window(xlim=xlim,ylim=c(0,1),asp=asp*xlim[2],xaxs="i",yaxs="i",bty="l")
	## draw axes
	axis(1,at=NULL,pos=0)
	axis(2,at=NULL,pos=0)
	#title(xlab=xlab)
	mtext(xlab, side=1, line=0,at=NA, adj=1)
	## draw the line
	lines(esf,type="s",col=col,lty=lty,lwd=lwd)

}


iwtp.esf2Plot <- function(esf1,esf2,xlim=NULL,xlab="",col=c("red","green"),lty=c(1,1),lwd=c(1,1),asp=0.6,...) {

	if(missing(xlim)||is.null(xlim))
		xlim <- c(0,max(esf1[,1],esf2[,1]))
	
	if(length(col)<2)
		col <- rep(col,times=2)
	if(length(lty)<2)
		lty <- rep(lty,times=2)
	if(length(lwd)<2)
		lwd <- rep(lwd,times=2)

	if(missing(xlab)||is.null(xlab))
		xlab = ""
	#if(missing(ylab)||is.null(ylab))
	#	ylab = ""	

	plot.new()
	plot.window(xlim=xlim,ylim=c(0,1),asp=asp*xlim[2],xaxs="i",yaxs="i",bty="l")
	## draw axes
	axis(1,at=NULL,pos=0)
	axis(2,at=NULL,pos=0)
	#title(xlab=xlab)
	mtext(xlab, side=1, line=0,at=NA, adj=1)
	lines(esf1,type="s",col=col,lty=lty,lwd=lwd)
	## draw the second line if not identical as first line
	if(!identical(esf1,esf2))
		lines(esf2,type="s",col=col[2],lty=lty[2],lwd=lwd[2],...)
	
}


## define generic function: esfPlot
esfPlot <- function(x,...) {UseMethod("esfPlot")}

plot.iwtp <- function(x,x2=NULL,xlim=NULL,xlab="",col="red",lty=1,lwd=1,asp=0.6,...) {
	if(missing(x2) || is.null(x2)) {
		iwtp.plot(x,xlim=xlim,xlab=xlab,col=col,lty=lty,lwd=lwd,asp=asp,...)
	}
	else
	{
		iwtp.plot2(x,x2,xlim=xlim,xlab=xlab,col=col,lty=lty,lwd=lwd,asp=asp,...)
	}
}

iwtp.plot <- function(sf1,xlim=NULL,xlab="",col=c("red","blue"),lty=c(1,1),lwd=c(1,1),asp=0.6,...) {
	sf.para <- sf1$sf$par
	if(class(sf1$sf)=="iwtp.interv.wb") {
		sf <- function(x) {sf.wb(x,sf.para[1],sf.para[2])}
		#cf <- mean.wb(sf1$sf$par)/mean.nonpar(sf1$esf$esf)
		cf <- mean.wb(sf1$sf$par)/sf1$esf$mean
	}
	
	if(class(sf1$sf)=="iwtp.interv.wemix") {
		sf <- function(x) {sf.wemix(x,sf.para[1],sf.para[2],sf.para[3],sf.para[4])}
		#cf <- mean.wemix(sf1$sf$par)/mean.nonpar(sf1$esf$esf)
		cf <- mean.wemix(sf1$sf$par)/sf1$esf$mean
	}
	
	if(length(col)<2)
		col <- rep(col,times=2)
	if(length(lty)<2)
		lty <- rep(lty,times=2)
	if(length(lwd)<2)
		lwd <- rep(lwd,times=2)
	
	if(is.null(xlim))
		x.value <- c(0:max(sf1$esf$esf[,1]))
	else
		x.value <- c(xlim[1]:xlim[2])
		
	if(missing(xlab)||is.null(xlab))
		xlab = ""
	#if(missing(ylab)||is.null(ylab))
	#	ylab = ""	

	xlim=range(x.value)
	plot.new()
	plot.window(xlim=xlim,ylim=c(0,1),asp=asp*xlim[2],xaxs="i",yaxs="i",bty="l")
	##draw axes
	axis(1,at=NULL,pos=0)
	axis(2,at=NULL,pos=0)
	#title(xlab=xlab)
	mtext(xlab, side=1,line=0,at=NA, adj=1)
	##draw empirical survival function
	lines(sf1$esf$esf,type="s",col=col,lty=lty,lwd=lwd)
		
	## draw survival function
	lines(sf(cf*x.value),col=col[2],lty=lty[2],lwd=lwd[2],...)	
}


iwtp.plot2 <- function(sf1,sf2,xlim=xlim,xlab="",col=c("red","blue","green","yellow"),lty=c(1,1,1,1),lwd=c(1,1,1,1),asp=0.6,...) {
	sf1.para <- sf1$sf$par
	sf2.para <- sf2$sf$par

	if(class(sf1$sf)=="iwtp.interv.wb") {
		sf.sf1 <- function(x) {sf.wb(x,sf1.para[1],sf1.para[2])}
		#cf.sf1 <- mean.wb(sf1$sf$par)/mean.nonpar(sf1$esf$esf)
		cf.sf1 <- mean.wb(sf1$sf$par)/sf1$esf$mean
	}
	
	if(class(sf1$sf)=="iwtp.interv.wemix") {
		sf.sf1 <- function(x) {sf.wemix(x,sf1.para[1],sf1.para[2],sf1.para[3],sf1.para[4])}
		#cf.sf1 <- mean.wemix(sf1$sf$par)/mean.nonpar(sf1$esf$esf)
		cf.sf1 <- mean.wemix(sf1$sf$par)/sf1$esf$mean
	}
	
	if(class(sf2$sf)=="iwtp.interv.wb") {
		sf.sf2 <- function(x) {sf.wb(x,sf2.para[1],sf2.para[2])}
		#cf.sf2 <- mean.wb(sf2$sf$par)/mean.nonpar(sf2$esf$esf)
		cf.sf2 <- mean.wb(sf2$sf$par)/sf2$esf$mean
	}
	
	if(class(sf2$sf)=="iwtp.interv.wemix") {
		sf.sf2 <- function(x) {sf.wemix(x,sf2.para[1],sf2.para[2],sf2.para[3],sf2.para[4])}
		#cf.sf2 <- mean.wemix(sf2$sf$par)/mean.nonpar(sf2$esf$esf)
		cf.sf2 <- mean.wemix(sf2$sf$par)/sf2$esf$mean
	}

	if(length(col)<4)
		col <- rep(col,times=4)
	if(length(lty)<4)
		lty <- rep(lty,times=4)
	if(length(lwd)<4)
		lwd <- rep(lwd,times=4)	
	
	if(is.null(xlim))
		x.value <- c(0:max(sf1$esf$esf[,1],sf2$esf$esf[,1]))
	else
		x.value <- c(xlim[1]:xlim[2])
		
	if(missing(xlab)||is.null(xlab))
		xlab = ""
	#if(missing(ylab)||is.null(ylab))
	#	ylab = ""	

	xlim=range(x.value)
	plot.new()
	plot.window(xlim=xlim,ylim=c(0,1),asp=asp*xlim[2],xaxs="i",yaxs="i",bty="l")
	## draw axes
	axis(1,at=NULL,pos=0)
	axis(2,at=NULL,pos=0)
	#title(xlab=xlab)
	mtext(xlab, side=1, line=0,at=NA, adj=1)
	lines(sf1$esf$esf,type="s",col=col,lty=lty,lwd=lwd)		

	##  draw survival function
	lines(sf.sf1(cf.sf1*x.value),col=col[2],lty=lty[2],lwd=lwd[2],...)

	## draw second empirical suvival funciton, if not identical as the first 
	if(!identical(sf1$esf,sf2$esf))	
		lines(sf2$esf$esf,type="s",xlim=xlim,col=col[3],lty=lty[3],lwd=lwd[3],...)
	## draw second survival function
	lines(sf.sf2(cf.sf2*x.value),col=col[4],lty=lty[4],lwd=lwd[4],...)
}



