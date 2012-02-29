


QQplotStepwise <- function(li,xlim=NULL,ylim=NULL,xlab="",ylab="",col="red",lty=1,lwd=1,asp=0.6,...) {
	## define local functions
	#erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
	cdfNr <- function(x) pnorm(x)
	
	n <- length(li)
	ml <- mean(li)
	l0 <- li - ml
	dy <- sort(l0)

	qt <- sapply(c(1:n/(n+1)),function(a) {f.root <- uniroot(function(x) {cdfNr(x)-a}, lower=min(dy), upper=max(dy)) ;
								return(f.root$root) } )

	if(missing(xlim)||is.null(xlim)) {
		rx <- range(qt)
		xlim=ceiling(abs(rx))*sign(rx) 
	}
	if(missing(ylim)||is.null(ylim)) {
		ry <- range(dy)
		ylim=ceiling(abs(ry))*sign(ry)
	}

	if(missing(xlab)||is.null(xlab))
		xlab = ""
	if(missing(ylab)||is.null(ylab))
	 	ylab = ""	
	
	## set graph window to square
	par(pty = "s")
	plot.new()
	asp=asp*diff(range(xlim))/diff(range(ylim))
	plot.window(xlim=xlim,ylim=ylim,asp=asp,xaxs="i",yaxs="i",bty="l")
	## get axes ticks
	xtiks <- pretty(xlim)
	ytiks <- pretty(ylim)

	## draw grid
	for( i in 1:length(xtiks))
		lines(rbind(c(xtiks[i],ytiks[1]),c(xtiks[i],ytiks[length(ytiks)])),col="lightgray",lwd=0.1,lty=1)

	for(i in 1:length(ytiks))
		abline(h=ytiks[i],col="lightgray",lwd=0.1,lty=1)	
	
	## draw axes
	axis(1,at=pretty(xlim),pos=ytiks[1])	
	axis(2,at=pretty(ylim),pos=xtiks[1])
	mtext(xlab, side=1, line=0,at=NA, adj=1)
	mtext(ylab, side=2, line=3,at=0, adj=1)
	#title(xlab=xlab,ylab=ylab)
	## draw the line
	lines(qt,dy,type="s",col=col[1],lty=lty[1],lwd=lwd[1],...)
	
	## add minor marks to x and y axis
	#axis(side = 1, at = seq(xlim[1], xlim[2], by = 1), labels = FALSE, tcl = -0.2) 
	#axis(side = 2, at = seq(ylim[1], ylim[2], by = 10), labels = FALSE, tcl = -0.2)
	return (invisible(cbind(qt,dy)))
}

##test
#n <- 200
#set.seed(5120); sd <- rndNr0(n)


QQplot2Stepwise <- function(l1,l2,xlim=NULL,ylim=NULL,xlab="",ylab="",col=c("red","blue"),lty=c(1,1),lwd=c(1,1),asp=0.6, ... ) {
	## define local functions
	#erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
	cdfNr <- function(x) pnorm(x)

	n1 <- length(l1) ; n2 <- length(l2)	
	m1 <- mean(l1)	; m2 <- mean(l2)
	ll1 <- l1 - m1    ; ll2 <- l2 - m2
	dy1 <- sort(ll1)   ; dy2 <- sort(ll2)

	qt1 <- sapply(c(1:n1/(n1+1)), function(a) {f.root <- uniroot(function(x) {cdfNr(x)-a}, lower=min(dy1), upper=max(dy1)) ;
								return(f.root$root)})
	qt2 <- sapply(c(1:n2/(n2+1)), function(a) {f.root <- uniroot(function(x) {cdfNr(x)-a}, lower=min(dy2), upper=max(dy2)) ;
								return(f.root$root)})

	if(length(col)<2)
		col <- rep(col,times=2)
	if(length(lty)<2)
		lty <- rep(lty,times=2)
	if(length(lwd)<2)
		lwd <- rep(lwd,times=2)

	if(missing(xlim)||is.null(xlim)) {
		rx <- range(qt1,qt2)
		xlim=ceiling(abs(rx))*sign(rx) 
	}
	if(missing(ylim)||is.null(ylim)) {
		ry <- range(dy1,dy2)
		ylim=ceiling(abs(ry))*sign(ry)
	}
	
	if(missing(xlab)||is.null(xlab))
		xlab = ""
	if(missing(ylab)||is.null(ylab))
		ylab = ""	

	## set graph window to square
	par(pty = "s")
	plot.new()
	asp=asp*diff(range(xlim))/diff(range(ylim))
	plot.window(xlim=xlim,ylim=ylim,asp=asp,xaxs="i",yaxs="i",bty="l")

	## get axes ticks
	xtiks <- pretty(xlim)
	ytiks <- pretty(ylim)

	## draw grid
	for( i in 1:length(xtiks))
		lines(rbind(c(xtiks[i],ytiks[1]),c(xtiks[i],ytiks[length(ytiks)])),col="lightgray",lwd=0.1,lty=1)

	for(i in 1:length(ytiks))
		abline(h=ytiks[i],col="lightgray",lwd=0.1,lty=1)	
	
	## draw axes
	axis(1,at=pretty(xlim),pos=ytiks[1])	
	axis(2,at=pretty(ylim),pos=xtiks[1])
	mtext(xlab, side=1, line=0,at=NA, adj=1)
	mtext(ylab, side=2, line=3,at=0, adj=1)
	#title(xlab=xlab,ylab=ylab)
	##draw the line
	lines(qt1,dy1,type="s",col=col[1],lty=lty[1],lwd=lwd[1],...)	

	## draw the line
	lines(qt2,dy2,type="s",col=col[2],lty=lty[2],lwd=lwd[2],...)
	
	## add minor marks to x and y axis
	##axis(side=1, at=seq(xlim[1],xlim[2],by=0.2), labels=FALSE,tcl=-0.2) 
	##axis(side=2, at=seq(ylim[1],ylim[2],by=2), labels=FALSE,tcl=-0.2)
}

##test
#n1 <- 100 ; n2 <- 100
#set.seed(5120); sd1 <- rndNr0(n1) ;
#set.seed(9760); sd2 <- rndNr0(n2) ;


## random real number generator
rndNr0 <- function(n) {
	m <- ceiling(n/2) ;
	l <- c() ;
	for(i in 1:m) {
		r <- runif(2) ;
		r1 <- sqrt(-2*log(r[1]))*cos(2*pi*r[2]) ;
		r2 <- sqrt(-2*log(r[1]))*sin(2*pi*r[2]) ;
		l <- c(l, r1, r2) ;
	}
	
	if(n != 2*m)
		l <- l[1:(2*m-1)] ;
	return(l)
}

## eof function rndNr0

