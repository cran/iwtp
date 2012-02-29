## function: esfl
## input: interval.R = T, based on the right end value of intervals and the corresponding ties
## output:
## note: empirical survival function
iwtp.esf <- function(data) {
	#if(interval.R) {
		order.wtp <- order(data[,1],data[,2])
		list.wtp <- cbind(data[,1][order.wtp],data[,2][order.wtp])
		esf.R <- listSF(list.wtp)
		out <- list(esf=esf.R, mean=mean.nonpar(esf.R))
		class(out) <- c("iwtp.esf")
		return(out)
	#}
}
## eof function esf


## function: mean.nonpar
## input:
## output:
## note: mean estimate for empirical survival function based on nonparametric estimate
## using higher prob
mean.nonpar <- function(sf) {
	n <- dim(sf)[1]
	m.wtp <- sf[1,1] + sum((sf[2:n,1]-sf[-n,1])*sf[-n,2])
	return(m.wtp)
}
## using lower prob
mean.nonpar.l <- function(sf) {
	n <- dim(sf)[1]
	m.wtp <- sf[1,1] + sum((sf[2:n,1]-sf[-n,1])*sf[-1,2])
	return(m.wtp)
}

## usign trape. 
mean.nonpar.t <- function(sf) {
	n <- dim(sf)[1]
	m.wtp <- sf[1,1] + sum((sf[2:n,1]-sf[-n,1])*(sf[-n,2]+sf[-1,2])/2)
	return(m.wtp)
}

calculateSFmean <- function(sfn, name, include.ends=F) {
	if(!include.ends) {
		sf.mean <- c(mean.nonpar.t(sfn))
		names(sf.mean) <- c(paste("mean.T.", name, sep=""))
	}
	else
	{
		sf.mean <- c(mean.nonpar(sfn), mean.nonpar.t(sfn),mean.nonpar.l(sfn) )
		names(sf.mean) <- c(paste("mean.R.","mean.T.","mean.L.", name, sep=""))
	}
	
	return(sf.mean)
	
}
## eof function mean.nonpar


## function: listSF 
## input:dati={{x1,t1},{x2,t2},...,{xn,tn}}
## output:list sf={{0,1},{x1,q1},...,{xn,0}}
## note:
listSF <- function(data) 
{
	if(is.null(data))
		stop("Specify at list `data' argument", call. = FALSE)
		
	if (!is.null(data) && is.matrix(data)) 
	{
		dt <- data; md <- length(dt[,1]); nd <- sum(dt[,2])
		df <- c(dt[,2]/nd)
		pr <- df 
		cp <- 1-cumsum(pr) 	
		mat <- matrix(dt[,1], byrow=TRUE)
		mat <- cbind(mat, cp)
    sf <- rbind(c(0,1), mat)
		colnames(sf) <- c("wtp","p")
		
		return (sf)
	}
}
## eof function listSF


