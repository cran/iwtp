# TODO: Add comment
# 
# Author: zhouw
###############################################################################


#The EM algorithm for interval censored data
# turnbull 1976
# this function are copied from package Icens (funciton EM)
turnbull<-function(A, pvec, maxiter = 500, tol = 1e-12)
{
	if( ncol(A)==2 && all(A[,2]>=A[,1]) ) {
		ml <- Maclist(A)
		intmap <- t(MLEintvl(A, ml))
		A <- Macmat(ml)$pmat
	}
	else
		intmap <- NULL
	i<-0
	notdone<-TRUE
	n<-ncol(A)
	Meps<-.Machine$double.eps
	if(missing(pvec))
		pvec <- apply(A, 1, sum)/sum(A)
	pvec<-rescaleP(pvec, Meps)
	while(i<maxiter && notdone) {
		i<-i+1
		dmat<-diag(pvec)
		t1<-dmat%*%A
		t2<-1/(t(A)%*%pvec)
		np<-rescaleP(as.vector(t1%*%t2)/n, Meps)
		if( sum(abs(np-pvec)) < tol )
			notdone<-FALSE
		pvec<-np
	}
	if( notdone )
		warning("EM may have failed to converge")
	ret <- list(pf=pvec, numiter=i,
			converge=!notdone, intmap=intmap)
	class(ret) <- "icsurv"
	return(ret)
}

auc.icsurv <- function(x) {
	if(class(x) != 'icsurv')
		stop("auc.icsurv only works with icsurv object")
	surv <- T
#	m<-dim(x$intmap)[[2]]
#	maxx<-1.1*max(x$intmap[!is.infinite(x$intmap)])
	effsupp<-x$intmap[,x$pf>0]
	effsupp<-cbind(c(0,0),effsupp)
	sigma <- x$sigma
	if( is.null(sigma) )
		sigma <- cumsum(x$pf)
	sigma<-c(0,sigma[x$pf>0])
	if (surv) sigma<-1-sigma
	
	delta <- diff(effsupp[2,])
	rectangles <- delta * sigma[-length(sigma)]
	mean <- sum(rectangles)
	mean
}

#auc.icsurv(tb)