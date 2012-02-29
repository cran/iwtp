## file: interv.R
## function: interv
## input: a data frame with self-selected intervals 
## output: a object of class "interv"
## usage: create divions based on the given set of self-selected intervals
## modified by zhouw, 2009-09-15
## modified by zhouw, 2010-03-23
## modified by zhouw, 2010-05-21
###################################################################################
interv <- function(data) {
	## bof local functions
	#
	# 	identify divisions
	##  return: matrix
	divisions <- function(x) {
		x <- x[,2:3]
		x <- sort(x)
		x <- union(x,c())  
		x <- rep(x,each=2)
		x <- x[2:(length(x)-1)]
		#print(x)
		return(matrix(x,ncol=2,byrow=TRUE))
	}
	
	#
	# 	division index for a given interval
	#   input x: a vector
	##  return vector
	divisionintv <- function(x,divs) {
		ndx.left <- which(divs[,1] == x[2])
		ndx.right <- which(divs[,2]== x[3])
		return (c(ndx.left, ndx.right))
	}
	
	#
	#	divison idex for all intervals
	#	input x: a matrix
	##  return: matrix
	ndx.divisions <- function(x, divs) {		
		mat <- t(apply(x,1,divisionintv,divs=divs))
		#print(mat)
		return (mat)
	}

	# 
	#	division size for all intervals
	##  return: matrix
	size.divisions <- function(ndx.divs){
		size <- apply(ndx.divs,1,function(x) {x[2]-x[1]+1})
		return (size)
	}
	
	#
	# aphj
	admissiblePairesHJ <- function(ndx.divs) {
		## ndx.divs: matrix, dim(ndx.divs)[2]=2
		################################
		
		nv <- max(ndx.divs[,2])
		paires <- c()
		for(j in 1:nv) {
			h <- as.vector(which(ndx.divs[,1]<=j&ndx.divs[,2]>=j))
			paires <- rbind(paires, cbind(h,j))
		}
		
		rownames(paires) <- c(1:dim(paires)[1])
		return(paires)
	}
		
	## eof local funcitons

	
	## bof function interv
	if(!is.matrix(data))
		data <- as.matrix(data)
	divs <- divisions(data)
	ndx <- ndx.divisions(data,divs)
	size.divs <- size.divisions(ndx)
	
	## allCh: list, a vector of division contained by a given interval
	allCh <- apply(ndx, 1, function(nd) {nd[1]:nd[2]}) 
	
	## aphj: all admissible  paires of h(stated interval) and j(division)
	aphj <- admissiblePairesHJ(ndx)
	
	## allDj, list, a vector of interval stated if the true value is located inside a given division
	nv <- dim(divs)[1]
	allDj <-sapply(1:nv, function(j){g<-aphj[which(aphj[,2]==j),1];names(g)=NULL;return(g)})

	## w0, vector, generting frequence of stated interval
	n.intv <- data[,4]
	w0 <- n.intv/sum(n.intv)
	
	whjBM <- calculateWhjBM(allDj, aphj, w0)
	
	vjR <- union(divs[,1], divs[,2])
	
	estprob.new <- estprobnew(data[,4],size=1000)
	
	out.list <- list(data=data,n.intv=n.intv,divs=divs,ndx.divs=ndx,size.divs=size.divs, vjR=vjR,
					aphj=aphj,allCh=allCh, allDj=allDj,whjBM=whjBM,w0=w0,estprob=estprob.new)
		
	class(out.list) <- "interv"
	return (out.list)
	
	## eof function interv
}

## define generic function Dofreq
Dofreq <- function(x,...) {
	UseMethod("Dofreq")
}

#
# S3 method
# estimates for frequency based on attractive
# argument x is a object of class 'divs.BM1', 
# argument intv is a object of class 'interv', an output from function interv 
## BM1: is the model of indifferent respondents with whj=1/dj
Dofreq.divs.BM1 <- function(x, intv, limits=list()) {
	# bof local functions
	tiesfordivisions <- function(intv) {
		mu <- length(intv$n.intv)
		kv <- dim(intv$divs)[1]
		ndx <- intv$ndx.divs
		pa <- cbind(intv$n.intv,ndx)
		mat <- t(apply(pa,1,function(xi,kl=0) {v <- rep(0,kl);v[xi[2]:xi[3]]<- 1;return(v)},kl=kv))
		rownames(mat) <- c(1:mu)
		#print(mat)
		return (mat)
	} 
	# eof local functions
	
	ties <- tiesfordivisions(intv) # 
	freqs <- t(t(ties)/colSums(ties))	
	out <- list(intv=intv,freq.mode="BM1",coef=list(),coef.limits=list(),freqs.divs=freqs)
	class(out) <- "iwtp.interv"
	return(out)
}


#
# S3 method
# estimates for frequency based on the last div
# argument x is a object of class 'divs.BM2', 
# argument intv is a object of class 'interv', an output from function interv 
## BM2 is the model of respondents who with vj containing their WTP-point 
## select to state uh in which vj is the last division interval.
## 
Dofreq.divs.BM2 <- function(x, intv,limits=list()) {
	# local functions
	tiesforHiDivision <- function(intv) {
		mu <- length(intv$n.intv)
		kv <- dim(intv$divs)[1]
		ndx <- intv$ndx.divs
		pa <- cbind(intv$n.intv,ndx)
		mat <- t(apply(pa,1,function(xi,kl=0) {v <- rep(0,kl);v[xi[3]]<- xi[1];return(v)},kl=kv))
		rownames(mat) <- c(1:mu)
		#print(mat)
		return (mat)
	} 
	# end of local functions
	
	ties <- tiesforHiDivision(intv)  # select only the last devision
	
	freqs <- t(t(ties)/colSums(ties))
	freqs[is.na(freqs)] <- 0
	
	out <- list(intv=intv,freq.mode="BM2",coef=list(),coef.limits=list(),freqs.divs=freqs)
	class(out) <- "iwtp.interv"
	return(out)
}

#
# S3 method
# argument x is a object of class 'divs', 
# argument intv is a object of class 'interv', an output from function interv 
## BM3: is the model of respondents who, with vj containing their wtp point, select
## to state uh, proportionally to achoring probabilities wh, h = 1, ...m
Dofreq.divs.BM3 <- function(x, intv, limits=list()) {
	# bof local functions
	tiesfordivisions <- function(intv) {
		mu <- length(intv$n.intv)
		kv <- dim(intv$divs)[1]
		ndx <- intv$ndx.divs
		pa <- cbind(intv$n.intv,ndx)
		mat <- t(apply(pa,1,function(xi,kl=0) {v <- rep(0,kl);v[xi[2]:xi[3]]<- xi[1];return(v)},kl=kv))
		rownames(mat) <- c(1:mu)
		#print(mat)
		return (mat)
	} 
	# eof local functions
	
	ties <- tiesfordivisions(intv) # 
	freqs <- t(t(ties)/colSums(ties))	
	out <- list(intv=intv,freq.mode="BM3",coef=list(),coef.limits=list(),freqs.divs=freqs)
	class(out) <- "iwtp.interv"
	return(out)
}

#
# S3 method
# argument x is an object of class 'divs.ranks', which is list containing the coefficient
# argument intv is an object of class 'interv', an output from function interv
# BM4: beta - type weighting of acnchoring probabilities wh, h = 1, ..., m
##
Dofreq.divs.BM4 <- function(x,intv, 
		limits=list(
				c1=c(lower=0.01,upper=100,init=0.1),
				c2=c(lower=0.01,upper=10,init=0.001))) {
	# bof local function
	ranksfordivisions <- function(intv) {
		mu <- length(intv$n.intv)
		kv <- dim(intv$divs)[1]
		ndx <- intv$ndx.divs
		ranks <- t(apply(ndx,1,function(xi,kl=0){v<-rep(0,kl);n<-xi[2]-xi[1]+1;v[xi[1]:xi[2]]<-1:n;return(v)},kl=kv))
		rownames(ranks) <- c(1:mu)
		#print(ranks)
		return (ranks)
	}
	
	do.frequencies <- function(th,ranks,coef.ranks=0,size.divs) {
		freq.div <- function(tr,th,size,coef=0) {
			# estimator of w
			w <- th/sum(th)
			# sum rank over division j, 
			#sum.tr <- sapply(size, function(x) sum(1:x))
			last.tr <- sapply(size, function(x) x)
			whj <- rep(0,length(th))
			
			# apply only to stated intervals u which contains division j
			idx <- which(tr>0)
			# calculate whj
			if(length(coef)==1)
				coef = c(coef,1)	
			
			#whj[idx] <- w[idx]*(1+coef[1]*(tr[idx]/sum.tr[idx])^coef[2])
			
			#
			rcj <- rep(0,length(th))
			#rcj[idx] <- tr[idx]/sum.tr[idx]
			rcj[idx] <- tr[idx]/last.tr[idx]
			
			whj[idx] <- rcj[idx]^(coef[1]-1) * ifelse((1-rcj[idx])==0,1,(1-rcj[idx])^(coef[2]-1))*w[idx]
			#
			
			sw <- sum(whj[idx])
			return(whj/sw) 
		}
		
		feq <- apply(ranks,2,freq.div,th=th,size=size.divs,coef=coef.ranks)
		return(feq)
	}
	
	# eof local functions
	
	ranks <- ranksfordivisions(intv)
	coef <- unlist(x)
	names(coef) <- names(x)
	freqs <- do.frequencies(intv$n.intv,ranks,coef.ranks=coef,intv$size.divs)
	
	out <- list(intv=intv,freq.mode="BM4",coef=coef,coef.limits=limits,freqs.divs=freqs)
	
	class(out) <- "iwtp.interv"
	return(out)		
}

# S3 method
# argument x is an object of class 'divs.ranks', which is a list containing the coefficient
# argument intv is an object of class 'interv', an output from function interv
# BM5: proportionally to wh*(1+c(rhj/rh.)))
##
Dofreq.divs.BM5 <- function(x,intv, 
		limits=list(
				c=c(lower=0.01,upper=100,init=0.1))) {
	# bof local function
	ranksfordivisions <- function(intv) {
		mu <- length(intv$n.intv)
		kv <- dim(intv$divs)[1]
		ndx <- intv$ndx.divs
		ranks <- t(apply(ndx,1,function(xi,kl=0){v<-rep(0,kl);n<-xi[2]-xi[1]+1;v[xi[1]:xi[2]]<-1:n;return(v)},kl=kv))
		rownames(ranks) <- c(1:mu)
		#print(ranks)
		return (ranks)
	}
	
	do.frequencies <- function(th,ranks,coef.ranks=0,size.divs) {
		freq.div <- function(tr,th,size,coef=0) {
			# estimator of w
			w <- th/sum(th)
			# sum rank over division j, 
			sum.tr <- sapply(size, function(x) sum(1:x) )
			whj <- rep(0,length(th))
			
			# apply only to stated intervals u which contains division j
			idx <- which(tr>0)
			# calculate whj
			
			# only one c
			whj[idx] <- w[idx]*(1+coef[1]*(tr[idx]/sum.tr[idx]))
			
			sw <- sum(whj[idx])
			return(whj/sw) 
		}
	
		feq <- apply(ranks,2,freq.div,th=th,size=size.divs,coef=coef.ranks)
		return(feq)
	}
	
	# eof local functions
	

	ranks <- ranksfordivisions(intv)
	coef <- unlist(x)
	names(coef) <- names(x)
	freqs <- do.frequencies(intv$n.intv,ranks,coef.ranks=coef,intv$size.divs)
	
	out <- list(intv=intv,freq.mode="BM5",coef=coef,coef.limits=limits,freqs.divs=freqs)
	
	class(out) <- "iwtp.interv"
	return(out)		
}


#
# S3 method
# estimate freqency based on ranks with two parameters. 
# argument x is an object of class 'divs.ranks', which is list containing the coefficient
# argument intv is an object of class 'interv', an output from function interv
##
Dofreq.divs.2 <- function(x,intv, 
		limits=list(
				c1=c(lower=0.01,upper=10000,init=0.1),
				c2=c(lower=0.01,upper=10,init=0.001))) {
	# bof local function
	ranksfordivisions <- function(intv) {
		mu <- length(intv$n.intv)
		kv <- dim(intv$divs)[1]
		ndx <- intv$ndx.divs
		ranks <- t(apply(ndx,1,function(xi,kl=0){v<-rep(0,kl);n<-xi[2]-xi[1]+1;v[xi[1]:xi[2]]<-1:n;return(v)},kl=kv))
		rownames(ranks) <- c(1:mu)
		#print(ranks)
		return (ranks)
	}
	
	do.frequencies <- function(th,ranks,coef.ranks=0,size.divs) {
		freq.div <- function(tr,th,size,coef=0) {
			# estimator of w
			w <- th/sum(th)
			# sum rank over division j, 
			sum.tr <- sapply(size, function(x) sum(1:x) )
			whj <- rep(0,length(th))
			
			# apply only to stated intervals u which contains division j
			idx <- which(tr>0)
			# calculate whj
			if(length(coef)==1)
				coef = c(coef,1)	
			
			whj[idx] <- w[idx]*(1+coef[1]*(tr[idx]/sum.tr[idx])^coef[2])
			
			#
			#rcj <- rep(0,length(th))
			#rcj[idx] <- tr[idx]/sum.tr[idx]
			#whj[idx] <- rcj[idx]^(coef[1]-1) * ifelse((1-rcj[idx])==0,1,(1-rcj[idx])^(coef[2]-1))*w[idx]
			#
			
			sw <- sum(whj[idx])
			return(whj/sw) 
		}
		
		feq <- apply(ranks,2,freq.div,th=th,size=size.divs,coef=coef.ranks)
		return(feq)
	}
	
	# eof local functions
	
	ranks <- ranksfordivisions(intv)
	coef <- unlist(x)
	names(coef) <- names(x)
	freqs <- do.frequencies(intv$n.intv,ranks,coef.ranks=coef,intv$size.divs)
	
	out <- list(intv=intv,freq.mode="2",coef=coef,coef.limits=limits,freqs.divs=freqs)
	
	class(out) <- "iwtp.interv"
	return(out)		
}



#
# S3 method
# estimate freqency based on the distance to midpoint of division. 
# argument x is an object of class 'divs.distances', which is a list containing the coefficients
# argument intv is an object of class 'interv', an output from function interv
##
Dofreq.divs.4 <- function(x,intv, 
		limits=list(
				c1=c(lower=0.01,upper=100,init=0.1),
				c2=c(lower=0.01,upper=10,init=0.001))) {
	
	#defin local function
	divisionsdistance <- function(intv) {
		ndx <- intv$ndx.divs
		divs <- intv$divs
		#size.divs <- intv$size.divs
		
		intv.range <- apply(ndx,1,function(x,divs) {divs[x[2],2]-divs[x[1],1]}, divs=divs)
		intv.left <- apply(ndx,1,function(x,divs) {divs[x[1],1]}, divs=divs)
		mid.divs <- apply(divs, 1, function(x) mean(x))
		
		mu <- dim(ndx)[1]
		kv <- dim(divs)[1]
		distanceToIntvLeft <- matrix(0,nrow=mu,ncol=kv)
		
		for(i in 1:mu) {
			distanceToIntvLeft[i, ndx[i,1]:ndx[i,2]] <- (mid.divs[ndx[i,1]:ndx[i,2]]-intv.left[i])/intv.range[i]
		}
		
		#rownames(ranks) <- c(1:mu)
		#print(distanceToIntvleft)
		return (distanceToIntvLeft)
	}
	
	do.frequencies <- function(th,distances,coef.ranks=0,size.divs) {
		# define local functions
		freq.div <- function(dist,th,size,coef=0) {
			# estimator of w
			w <- th/sum(th)
			# sum rank over division j, 
			#sum.tr <- sapply(size, function(x) sum(1:x) )
			whj <- rep(0,length(th))
			
			# apply only to stated intervals u which contains division j
			idx <- which(dist>0)
			# calculate whj
			if(length(coef)==1)
				coef = c(coef,1)	
			
			#whj[idx] <- w[idx]*(1+coef[1]*(tr[idx]/sum.tr[idx])^coef[2])
			
			#
			xcj <- rep(0,length(th))
			xcj[idx] <- dist[idx]
			whj[idx] <- xcj[idx]^(coef[1]-1) * ifelse((1-xcj[idx])==0,1,(1-xcj[idx])^(coef[2]-1))*w[idx]
			#
			
			sw <- sum(whj[idx])
			return(whj/sw) 
		}
		# end local functions
		
		feq <- apply(distances,2,freq.div,th=th,size=size.divs,coef=coef.ranks)
		return(feq)
	}
	
	# eof local functions
	
	#ranks <- ranksfordivisions(intv)
	#freqs <- do.frequencies(intv$n.intv,ranks,coef.ranks=x$coef,intv$size.divs)
	
	distances <- divisionsdistance(intv)
#		freqs <- calFrequencies(intv$n.intv,distances,coef.ranks=x$coef,intv$size.divs)
	coef <- unlist(x)
	names(coef) <- names(x)
	freqs <- do.frequencies(intv$n.intv,distances,coef.ranks=coef,intv$size.divs)
	out <- list(intv=intv,freq.mode="4",coef=coef,coef.limits=limits,freqs.divs=freqs)
	
	class(out) <- "iwtp.interv"
	return(out)		
}


## 
calculateWhjBM <- function(allDj,aphj,w0){
	## allDj: list
	## aphj: matrix, dim(aphj)[2]=2
	## w0: vector
	## BM:3
	#######################################
	whj <- apply(aphj,1, function(hj) {
				h <- hj[1]
				j <- hj[2]
				sDj <- allDj[[j]]
				tWj <- w0[sDj]
				sj <- sum(tWj)
				w <- ifelse(sj>0, w0[h]/sj, 0)
				return(c(h,j,w))
			})
	rownames(whj) <- c("h","j","whj")
	return(t(whj))
	
}


## eof function interv

print.iwtp.interv <- function(x, ...) {

	x <- x$intv
    cat("Respondent no.: ", sum(x$n.intv), "\n")
	cat(sprintf("Self-selected interval no.: %d\n",length(x$n.intv)))

	cat("\n")
	cat(sprintf("Division intervals: %d in total\n", dim(x$divs)[1]))
	n <- dim(x$divs)[1]
	count <- 0
	for(i in 1:n){
		if(i != n)
			cat(sprintf("(%d, %d], ",x$divs[i,1],x$divs[i,2]))
	    else
		  cat(sprintf("(%d, %d]",x$divs[i,1],x$divs[i,2]))
	  
		count <- count+1
		if(count>4) {
			cat("\n")
			count<-0
		}
	}

	cat("\n")
	#cat("Frequences of division\n")
	#options(digits=4)
	#for(i in 1:dim(x$freq.divs)[1]){
	#	cat(i, ": ")
	#	cat(subset(x$freq.divs[i,], x$freq.divs[i,]>0))
	#	cat("\n")
	#}

	#cat("\n\n")

}


defaultPara <- function(mode="1"){
	
	freq.fun <- paste("Dofreq","divs", mode, sep=".")
	fm <- formals(freq.fun)
	fm$limits
}

#test
#duhT <- c(c(1,4,11,4), c(2,1,7,5),	c(3,6,14,2), c(4,3,9,1))
#duhT <- matrix(duhT, ncol=4, byrow=TRUE)
#iwtp1<-interv(data.frame(duhT))

## plots

plot.interv <- function(x,plot.type="itv",xlab="",ylab="",...) {
	
	if(!tolower(plot.type) %in% c("itv","phj"))
		stop("plot.tyhpe must be itv, phj !", call. = FALSE)		
	if(plot.type=='itv') {
		interv.intv.plot(x,xlab="",ylab="",...)
	}
	else {
		interv.phj.plot(x,xlab="",ylab="",...)
	}
}


interv.intv.plot <- function(intv,xlab="",ylab="",title="",...) {
	data <- intv$data
	plot.new()
	x.lim <- c(0,max(data[,3]))
	y.lim <- c(0,dim(data)[1]+5)
	plot.window(xlim=x.lim,ylim=y.lim,xaxs='i',yaxs="i",bty="l")
	axis(1,at=pretty(x.lim))
	axis(2,at=pretty(y.lim))
	#plot(data[,2:3])
	for(i in 1:dim(data)[1]) {
		#print(c(data[i,3],i))
		lines(c(data[i,2],c(data[i,3])),c(i,i),col='red',lwd=data[i,4]/5)
		#print(c(data[i,2],i),c(data[i,3],i))
	} 	
	title(main=title)
	mtext(xlab,side=1,line=2,at=NA,adj=1)
	mtext(ylab,side=2,line=2,at=NA,adj=1)
}

interv.phj.plot <- function(intv,xlab="",ylab="",title="",...) {
	#intv
	aphj <- intv$aphj
	data <- intv$data
	x.lim <- c(0,dim(data)[1]+5)
	y.lim <- c(0,dim(intv$divs)[1]+2)
	plot(aphj,type='p',pch=20,lwd=1,bty='l',xaxs='i',yaxs="i",xlim=x.lim,ylim=y.lim,xlab="",ylab="")
	title(main=title)
	mtext(xlab,side=1,line=2,at=NA,adj=1)
	mtext(ylab,side=2,line=2,at=NA,adj=1)
}
