# S-plus/R functions to determine the NPMLE of the distribution for interval-
# censored event time.
# Copyright 1998-2000 Alain Vandal and Robert Gentleman
# University of Auckland
# These functions should work, but their intent is to illustrate the
# concepts involved.  Functions are provided as is, with no guarantee.
# Redistribute freely, without undocumented modification & without charge.
# Queries to vandal@stat.auckland.ac.nz or rgentlem@stat.auckland.ac.nz.

# Returns a list of maximal antichains from a list of real valued intervals
# Arguments:  intvls:  n x 2 matrix;first column contains left endpoints
#				second column contains right endpoints
# Returned value:   list of length m (magnitude of underlying interval order)
#		      - each entry corresponds to one maximal antichain
#		      - each entry contains the row numbers of all intervals
#					  belonging to the maximal antichains
#		      - maximal antichains occur in the list in their natural
#					  linear ordering
# Known bugs: In R, will issue some ignorable warnings if there is
#right-censored data (with an "Inf" right endpoint).
Maclist <- function(intvls, Lopen=TRUE, Ropen=FALSE)
{
    m <- dim(intvls)[1]
    id <- 1:m
    or <- order(intvls[,1])
    maclist <- NULL
    curmac <- id[or[1]]
    minend <- intvls[curmac,2]
    for (i in 2:m) {
    	curintvl <- id[or[i]]
        if( intvls[curintvl,1]>minend ||
           ((Lopen || Ropen) &&  intvls[curintvl,1]==minend ) ) {
                                        # New maximal antichain
            maclist <- c(maclist,list(curmac))
            oldmac <- curmac
            curmac <- NULL
            for (j in 1:length(oldmac))
                if ( intvls[curintvl,1]<intvls[oldmac[j],2] ||
                 (!Lopen && !Ropen &&
                  intvls[curintvl,1]==intvls[oldmac[j],2]) )
                    curmac <- c(curmac,oldmac[j])
            curmac <- c(curmac,curintvl)
            minend <- min(intvls[curmac,2])
        } else {
            curmac <- c(curmac,curintvl)
            minend <- min(minend,intvls[curintvl,2]) }
    }
    c(maclist,list(curmac))
}

# Returns the clique matrix and Petrie pairs of an interval order
# given its list of maximal antichains Arguments: ml: list of maximal
# antichains as returned by Maclist Returned value: object containing
# # - pmat: clique matrix of the underlying interval order, # rows are
# ordered according to the linear ordering of # the maximal antichains
# # - ppairs: Petrie pairs indicate the first and last # maximal
# antichains to which each elements belongs

Macmat <- function(ml)
{
    temp <- NULL
    m <- length(ml)
    for (i in 1:m)
        temp <- c(temp,ml[[i]])
    temp <- sort(unique(temp))
    n <- length(temp)
    ppairs <- matrix(0,2,n)
    retmat <- matrix(0,m,n)
    for (i in 1:m) {
        for (j in ml[[i]]) {
            if (ppairs[1, j]==0)
                ppairs[1, j] <- i
            ppairs[2, j] <- i
        }
        retmat[i, ml[[i]]] <- 1
    }
    dimnames(ppairs) <- list(c("Start","End"),temp)
    dimnames(retmat) <- list(NULL,temp)
    ret <- list(pmat = retmat, ppairs = ppairs)
    class(ret) <- "petrie"
    return(ret)
}

# Produce the mapping of the maximal antichains to their real interval
# representation for an interval order given by real-valued intervals.
# Arguments:	intvls:	see Maclist
#		ml:	list of maximal antichains for the intervals as
#			returned by Maclist
# Returned values:  matrix m x 2 containing the mapping row-wise
#		(1rst row corresponds to 1rst maximal antichains, etc.)
#		     m is the number of maximal antichains,
#		     1rst column contains left endpoints of the mapping
#		     2nd column contains right endpoints of the mapping

MLEintvl <- function(intvls, ml=Maclist(intvls))
{
    if( ncol(intvls) != 2 || any(intvls[,2] < intvls[,1]) )
        stop("only one dimensional intervals can be handled")
    m <- length(ml)
    ret <- matrix(0, m, 2)
    for(i in 1:m) {
        LL <- min(intvls)
        RR <- max(intvls)
        for(j in 1:length(ml[[i]])) {
            LL <- max(LL, intvls[ml[[i]][j], 1])
            RR <- min(RR, intvls[ml[[i]][j], 2])
        }
        ret[i,  ] <- c(LL, RR)
    }
    ret
}

# Pool monotone groups algorithm
#  Adapted from Y.L. Zhang & M.A. Newton (1997)
# (http://www.stat.wisc.edu/~newton/newton.html)
# Isotonizes a weighted and ordered set of values
# Arguments:	est:	the list of values
#				ww:		their weights
# Returned values: object containing
#		- est:  isotonized estimates
#		- ww:	weights of the isotonized estimates
#		- poolnum:  number of values pooled in the current
#					  estimate
#		- passes:  number of passes which were required to
#					  isotonize the list

PMGA<-function(est,ww=rep(1,length(est)))
{
    curm<-length(est)
    poolnum<-rep(1,curm)
    passes<-0
    iso<-FALSE
    while (!iso) {
        iso<-TRUE
    	poolind<-1
    	curind<-1
    	while (curind<=curm) {
            groupstart<-curind
            while (curind<curm && est[curind+1]<est[curind])
                curind<-curind+1
            iso<-poolind == curind
            est[poolind]<-sum(ww[groupstart:curind]*est[groupstart:curind])/sum(ww[groupstart:curind])
            ww[poolind]<-sum(ww[groupstart:curind])
            poolnum[poolind]<-sum(poolnum[groupstart:curind])
            poolind<-poolind+1
            curind<-curind+1
        }
        curm<-poolind-1
        passes<-passes+1
    }
    return(list(est = est[1:curm], ww = ww[1:curm], poolnum = poolnum[1:curm],
           passes = passes))
}

# Returns the (unsmoothed) NPMLE of the distribution function on the maximal
# antichains of interval censored survival data.
# The algorithm is adapted from Wellner & Zhan (1997).
# Arguments:
#	- A:  clique matrix of the data (only necessary argument)
#	- EMstep:  boolean determining whether an EM-step will be taken
#	  at each iteration
#	- ICMstep: boolean determining whether an ICM step will be taken
#	  at each iteration
#	- checkbnds: make sure that isotonization step does not wash out
#	  essential maximal antichains by using self-consistent bounds
#	- keepiter: boolean determining whether to keep the iteration
#	  states
#	- eps: maximal L1 distance between successive estimates before
#	  stopping iteration
#	- maxiter:  maximal number of iterations to perform before
#	  stopping
# Returned values:  object containing
#	- sigma:  NPMLE of the survival function on the maximal
#	  antichains
#	- weights:  diagonal of the likelihood function's second
#	  derivative
#	- lastchange:  vector of differences between the last two
#	  iterations
#	- numiter:  total number of iterations performed
#	- iter: (only present if keepiter is true) states of sigma during
#	  the iteration

EMICMmac <- function(A, EMstep=TRUE, ICMstep=TRUE, keepiter=FALSE, tol=1e-7,
                     tolbis=1e-7,maxiter=1000)
{
    if (!EMstep && !ICMstep) {
        print("One of EMstep or ICMstep must be true.")
        return(NULL)
    }
    Meps<-.Machine$double.eps
    m<-dim(A)[1]
    n<-dim(A)[2]
    tA<-t(A)
    if (m==1) {
        ret<-NULL
        ret$sigma<-1
        ret$weights<-n
        ret$lastchange<-0
        ret$numiter<-0
        return(ret)
    }

    WW<-matrix(0,m,n)
    for (i in 1:(m-1))
        WW[i,]<-A[i,]-A[i+1,]
    WW[m,]<-A[m,]
    sigma<-cumsum(apply(A,1,sum)/sum(A))
    numiter<-0
    oldsigma<-rep(-1,m)
    if (keepiter) iter<-sigma
    while (max(abs(oldsigma-sigma))>tol && numiter<=maxiter) {
        oldsigma<-sigma
        if (EMstep) {
            pvec<-diff(c(0,sigma))
            temp<-sweep(A,1,pvec,FUN="*")
            if (sum(apply(temp,2,sum)==0)==0) {
                pvec<-apply(sweep(temp,2,apply(temp,2,sum),
                                  FUN="/"),1,sum)
                sigma<-cumsum(pvec)/sum(pvec)
            }
            if (keepiter) iter<-rbind(iter,sigma)
        }
        if (ICMstep) {
            Wps<-1/(t(WW)%*%sigma)
            weights<-(abs(WW)%*%Wps^2)
            increment<-as.vector((WW%*%Wps)/weights)
            sigma<-sigma+increment
            sigma[m]<-1
            if (keepiter) iter<-rbind(iter,sigma)
            temp<-PMGA(sigma[-m],weights[-m])
            poolnum<-c(0,cumsum(temp$poolnum))
            for (i in 2:length(poolnum))
                for (j in (poolnum[i-1]+1):poolnum[i])
                    sigma[j]<-temp$est[i-1]
            if (keepiter) iter<-rbind(iter,sigma)
           # Implementing Jongbloed's correction through bisection
            temp<-c(0,sigma)
            pvec<-diff(c(0,oldsigma))
            ndir<-diff(c(0,temp[2:(m+1)]))-pvec
            pvec<-Bisect(tA,pvec,ndir,Meps,tolbis=1e-7)
            sigma<-cumsum(pvec)
            if (keepiter) iter<-rbind(iter,sigma)
        }
        numiter<-numiter+1
    }
    if (numiter == maxiter)
        warning("EM/ICM may have failed to converge.")
    pf<-diff(c(0,sigma))
    ret<-list(sigma=sigma,pf=pf,llk=sum(log(t(A)%*%pf)),
              weights=as.vector(weights),lastchange=sigma-oldsigma,
              numiter=numiter,eps=tol)
    if (keepiter) {
        if (EMstep && ICMstep)
            dimnames(iter)<-list(c("Seed",rep(c("EM","Fisher","PMGA","Bisect"),
                                              numiter)),NULL)
        else if (EMstep)
            dimnames(iter)<-list(rep("EM",numiter+1),NULL)
        else
            dimnames(iter)<-list(c("Seed",rep(c("Fisher","PMGA"),
                                              numiter)),NULL)
        ret$iter<-iter
    }
    ret
}


# Returns the distribution function NPMLE for interval censored data, with
# information regarding its real-line support
# Arguments:  - intvls: list of intervals as per Maclist
#			  - all other arguments:  see EMICMmac
# Returned values: object containing
#	      - all information as per returned value of EMICMmac
#	      - pf:  probability function on the maximal antichains
#	      - intmap:  real interval mapping for the mass of the
#		    NPMLE;  the Groeneboom-Wellner estimate is derived
#		    by assigning all the mass of the NPMLE on the maximal
#	            antichain to the right endpoint of the corresponding
#                   interval.
#	      - class:  value "icsurv"

EMICM <- function(A, EMstep=TRUE, ICMstep=TRUE, keepiter=FALSE, tol=1e-7,
                  maxiter=1000)
{
    if( ncol(A) == 2 && all(A[,2]>=A[,1]) ) {
        ml<-Maclist(A)
        intmap <- t(MLEintvl(A, ml))
        A <- Macmat(ml)$pmat
    }
    else
        intmap <- NULL
    temp<-EMICMmac(A, EMstep=EMstep, ICMstep=ICMstep,
                   keepiter=keepiter,tol=tol,maxiter=maxiter)
    if (is.null(temp)) return(NULL)
    class(temp) <- "icsurv"
    temp$intmap <- intmap
    temp
}

# Plotting an icsurv class object.
# Arguments:
#			- x:  an object returned by EMICM
#			- type:		"eq" equivalence class
#						"gw" Groeneboom-Wellner estimate
#						"lc" left-continuous estimate
# Returned value:  none;  a plot is produced on the current graphics device.
plot.icsurv <- function(x, type="eq", surv=FALSE, bounds=FALSE, shade=3,
                        density=30, angle=45, lty=1, new=TRUE, xlab="Time",
                        ylab="Probability", main="GMLE",ltybnds=2,...)
{
    if(!inherits(x, "icsurv" ) ) {
	stop("plot.icsurv only works for icsurv objects")
    }
    if( is.null(x$intmap) )
        stop("the object does not contain the real representation")

    m<-dim(x$intmap)[[2]]
    maxx<-1.1*max(x$intmap[!is.infinite(x$intmap)])
    effsupp<-x$intmap[,x$pf>0]
    effsupp<-cbind(c(0,0),effsupp)
    sigma <- x$sigma
    if( is.null(sigma) )
        sigma <- cumsum(x$pf)
    sigma<-c(0,sigma[x$pf>0])
    if (surv) sigma<-1-sigma
    if (new) {
        if (missing(main)) {
            if (surv) {
                if (type=="gw")
                    main<-"Survival GWMLE function"
                else if (type=="lc")
                    main<-"Left-continuous survival GMLE function"
                else
                    main<-"GMLE survival equivalence class"
            } else {
                if (type=="gw")
                    main<-"Distribution GWMLE function"
                else if (type=="lc")
                    main<-"Left-continuous distribution GMLE function"
                else
                    main<-"GMLE distribution equivalence class"
            }
        }
        plot(c(0,maxx),c(0,1),xlab=xlab,ylab=ylab,main=main,type="n",...)
    }
    effm<-length(sigma)-1
    if((type=="gw" && !surv) || (type=="lc" && surv)) {
        for (i in 1:effm) {
            segments(effsupp[2,i],sigma[i],effsupp[2,i+1],sigma[i],lty=lty)
            segments(effsupp[2,i+1],sigma[i],effsupp[2,i+1],sigma[i+1],lty=lty)
        }
        if (sum(is.infinite(effsupp))==0)
            segments(effsupp[2,effm+1],1,maxx,1)
    }
    else if ((type=="lc" && !surv) ||( type=="gw" && surv)) {
        for (i in 1:effm) {
            segments(effsupp[1,i],sigma[i],effsupp[1,i+1],sigma[i],lty=lty)
            segments(effsupp[1,i+1],sigma[i],effsupp[1,i+1],sigma[i+1],lty=lty)
        }
        if (sum(is.infinite(effsupp))==0)
            segments(effsupp[1,effm+1],1,maxx,1)
    }

    else {
        if (shade==1) {
            for (i in 1:effm) {
                    segments(effsupp[2,i],sigma[i],effsupp[1,i+1],sigma[i],lty=lty)
                }
        } else if (shade==2) {
            for (i in 1:effm) {
                segments(effsupp[2,i],sigma[i],effsupp[1,i+1],sigma[i],lty=lty)
                    polygon(c(effsupp[2,i+1],effsupp[2,i+1],effsupp[1,i+1],effsupp[1,i+1]),
                            c(sigma[i],sigma[i+1],sigma[i+1],sigma[i]),
                            border=TRUE,density=0,lty=2)
            }
        } else {
            for (i in 1:effm) {
                segments(effsupp[2,i],sigma[i],effsupp[1,i+1],sigma[i],lty=lty)
                                        #dirty check on R vs. S, since density not applied in R
                if (is.null(version$language))
                    polygon(c(effsupp[2,i+1],effsupp[2,i+1],effsupp[1,i+1],effsupp[1,i+1]),
                            c(sigma[i],sigma[i+1],sigma[i+1],sigma[i]),
                            border=FALSE,density=density,angle=angle)
                else
                    polygon(c(effsupp[2,i+1],effsupp[2,i+1],effsupp[1,i+1],effsupp[1,i+1]),
                            c(sigma[i],sigma[i+1],sigma[i+1],sigma[i]),
                            border=FALSE,col="green")
            }
                if (bounds) {
                    bndhi<-cbind(c(0,x$intmap[1,]),c(x$intmap[1,],maxx))
                    bndlo<-cbind(c(0,x$intmap[2,]),c(x$intmap[2,],maxx))
                    if (surv) bndval<-1-bndval
                    segments(bndlo[,1],x$bounds[1,],bndlo[,2],x$bounds[1,],lty=ltybnds)
                    segments(x$intmap[2,],x$bounds[1,-(m+1)],x$intmap[2,],x$bounds[1,-1],lty=ltybnds)
                    segments(bndhi[,1],x$bounds[2,],bndhi[,2],x$bounds[2,],lty=ltybnds)
                    segments(x$intmap[1,],x$bounds[2,-(m+1)],x$intmap[1,],x$bounds[2,-1],lty=ltybnds)
                }
        }
        if (sum(is.infinite(effsupp))==0)
            if (!surv)
                segments(effsupp[2,effm+1],1,maxx,1,lty=lty)
            else
                segments(effsupp[2,effm+1],0,maxx,0,lty=lty)
    }
}


####################################
# MIXTURE METHODS for interval censored data survival estimation
# VEM, ISDM, PGM
# All require argument A, the clique matrix of the data, so that a typical call would be
# VEM(Macmat(Maclist(brcm))$pmat)
##################################
Bisect <- function(tA, pvec, ndir, Meps, tolbis=1e-7)
{
    etainv<-1/(tA%*%pvec)
    bot<-0
    top<-1
    mult<-tA%*%ndir
    dbot<-sum(etainv*mult)
    ptop<-rescaleP(pvec+top*ndir, Meps)
    pbot<-pvec
    dtop<-sum(mult/(tA%*%ptop))
    done<-FALSE
    while( !done ) {
        if( sign(dbot)*sign(dtop) > 0 || top-bot<tolbis ) {
            ltop<-sum(log(tA%*%ptop))
            lbot<-sum(log(tA%*%pbot))
            if( lbot > ltop )
                pnew<-rescaleP(pvec+bot*ndir, Meps)
            else
                pnew<-rescaleP(pvec+top*ndir, Meps)
            done<-TRUE
        }
        else {
            mid<-(bot+top)/2
            pmid<-rescaleP(pvec+mid*ndir, Meps)
            dmid<-sum(mult/(tA%*%pmid))
            if( dmid*dtop < 0 ) {
                bot<-mid
                dbot<-dmid
                pbot<-pmid
            }
            else {
                top<-mid
                dtop<-dmid
                ptop<-pmid
            }
        }
    }
    pnew
}

##################################
# Bohning's Vertex Exchange Method
##################################
VEM <- function(A, pvec, maxiter = 500, tol = 1e-7, tolbis = 1e-7,
                keepiter=FALSE)
{
    if( ncol(A) == 2 && all(A[,2]-A[,1] >= 0) ) {
        ml <- Maclist(A)
        intmap <- t(MLEintvl(A, ml))
        A <- Macmat(ml)$pmat
    }
    else
        intmap <- NULL
    m<-dim(A)[1]
    n<-dim(A)[2]
    tA<-t(A)
    Meps<-.Machine$double.eps
    if(missing(pvec))
        pvec <- apply(A, 1, sum)/sum(A)
    pvec<-rescaleP(pvec, Meps)
    if (keepiter)
        iter<-pvec
    m<-length(pvec)
    Linv<-1/(tA%*%pvec)
    lnew<--sum(log(Linv))
    finished<-FALSE
    i<-0
    while( i<maxiter && !finished ) {
        i<-i+1
        dvec<-A%*%Linv
        mind<-min(dvec[pvec>0])
        minj<-match(mind,dvec)
        maxd<-max(dvec)
        maxj<-match(maxd,dvec)
        ndir<-rep(0,m)
        ndir[minj]<- -1
        ndir[maxj]<-1
        ndir<-ndir*pvec[minj]
        pold<-pvec
        pvec<-Bisect(tA,pvec,ndir,Meps,tolbis=tolbis)
        lold<-lnew
        Linv<-1/(tA%*%pvec)
        lnew<--sum(log(Linv))
        if (keepiter)
            iter<-cbind(iter,pvec)
        finished<-(sum(abs(pvec-pold))<tol && abs(lnew-lold)<tol )
    }
    if( !finished )
        warning("VEM may have failed to converge")
    ret<-list(pf=pvec, lval=lnew, numiter=i, intmap=intmap,
              converge=finished )
    if (keepiter)
        ret$iter<-iter
    class(ret) <- "icsurv"
    ret
}

#################################
# Lesperance & Kalbfleisch ISDM
#################################

ISDM <- function(A, pvec, maxiter = 500, tol = 1e-07, tolbis = 1e-08,
               verbose = FALSE )
{
    if( ncol(A) == 2 && all( A[,2]-A[,1] >= 0 ) ) {
	ml <- Maclist(A)
        intmap <- t(MLEintvl(A, ml))
	A <- Macmat(ml)$pmat
    }
    else
        intmap <- NULL
    n <- dim(A)[2]
    m <- dim(A)[1]
    Meps<-.Machine$double.eps
    tA<-t(A)
    if(missing(pvec))
        pvec <- apply(A, 1, sum)/sum(A)
    pvec<-rescaleP(pvec, Meps)
    LL <- as.vector(tA %*% pvec)
    dd <- A %*% (1/LL) - n
    ind <- dd > 0
    numiter <- 0
    finished <- FALSE
    eps0 <- 0.5  #or some other number
    while(!finished && (numiter < maxiter) ) {
        numiter <- numiter + 1
        k <- sum(ind)
	eps0 <- 1-k/m  #or some other number
        aug <- rbind(LL, A[ind, ])
        taug <- t(aug)
        p1 <- c(eps0, (1-eps0)*(pvec[ind]/sum(pvec[ind])))
        rval <- VEM(aug, p1, tol=1e-5/numiter, tolbis=1e-5/numiter)
        epsilon <- rval$pf
        curlike <- prod(LL)
        LL <- as.vector(taug%*% epsilon)
        if(verbose) {
            print(paste("iteration: ",numiter))
            print(paste("VEM iterations: ", rval$numiter))
            print(epsilon)
            print(round(pvec,4))
            print(sum(log(LL)))
        }
        pold <- pvec
        pvec <- as.vector(solve(A %*% tA, A %*% LL))
        pvec <- rescaleP(pvec, Meps)
        finished <- sum(abs(pvec - pold)) < tol
        dd <- A %*% (1/LL) - n
        ind <- dd >= 0
    }
    if( !finished )
        warning("ISDM may have failed to converge")
    ret <- list(pf=pvec, numiter=numiter, intmap=intmap,
                converge=finished)
    class(ret) <- "icsurv"
    return(ret)
}


#############
#rescaleP is a function that rescales a prob vector so that elements
# that are negative or less than machine epsilon are set to zero.
###########
rescaleP <- function(pvec, tiny)
{
    pvec<-ifelse(pvec<tiny,0,pvec)
    pvec<-pvec/sum(pvec)
    return(pvec)
}

########################
# Projected gradient method
# Search direction is based on all directional derivatives on the
# caveat that support points with 0 mass and negative directional
# derivatives are excluded from the direction.
# Alain: told does not seem to be used???
########################
PGM <- function(A, pvec, maxiter = 500, tol = 1e-07, told=2e-5,
                tolbis = 1e-08, keepiter=FALSE)
{
    if( ncol(A) == 2 && all(A[,2]>=A[,1]) ) {
        ml <- Maclist(A)
        intmap <- t(MLEintvl(A, ml))
        A <- Macmat(ml)$pmat
    }
    else
        intmap <- NULL
    m<-dim(A)[1]
    n<-dim(A)[2]
    tA<-t(A)
    Meps<-.Machine$double.eps
    if(missing(pvec))
        pvec <- apply(A, 1, sum)/sum(A)
    Linv<-1/(tA%*%pvec)
    lnew<-sum(-log(Linv))
    dd <- A %*% Linv
    if (keepiter)
        iter<-pvec
    i <- 0
    finished <- FALSE
    while (i<maxiter && !finished) {
        i <- i + 1
        pvec<-rescaleP(pvec, Meps)
        # Limit directional derivative vector to possible
        # directions of motion (!zero)
        zero<-dd<n & pvec<=0
        # Get projected gradient
        # Lam can be replaced with another pd matrix to get
        # another member of Wu's class of search directions
        dlam<-dd[!zero]
        Lam<-diag(rep(1,m-sum(zero)))
        dir<-rep(0,m)
        dir[!zero]<-(m-sum(zero))*dlam-sum(dlam)
        # Renormalize for numerical stability
        dir<-dir/max(abs(dir))
        # Obtain maximum distance from current vector
        # in direction of dir
        pos<-dir>0
        neg<-dir<0
        pma<-min(c(-pvec[neg]/dir[neg],(1-pvec[pos])/dir[pos]))
        dir<-pma*dir
        # Setup the bisection
        pold<-pvec
        pvec<-Bisect(tA,pvec,dir,Meps,tolbis=tolbis)
        if (keepiter)
            iter<-cbind(iter,pvec)
        lold<-lnew
        Linv<-1/(tA%*%pvec)
        lnew<-sum(-log(Linv))
        dd <- A %*% (1/(tA%*%pvec))
        finished <- (sum(abs(pvec - pold)) < tol && abs(lnew-lold)<tol)
    }
    if ( !finished )
        warning("PGM may have failed to converge.")
    eps<-c(tol,told,tolbis)
    names(eps)<-c("tol","told","tolbis")
    ret<-list(pf=pvec,sigma=cumsum(pvec),lval=lnew,clmat=A,
              method="MPGM", lastchange=pvec-pold, numiter=i,
              eps=eps, converge=finished, intmap=intmap)
    if (keepiter)
        ret$iter<-iter
    class(ret) <- "icsurv"
    ret
}

# Returns the (unsmoothed) NPMLE of the distribution function on the maximal
# antichains of interval censored survival data.
# The algorithm is adapted from Wellner & Zhan (1997).
# Arguments:
#	- A:  clique matrix of the data (only necessary argument)
#	- EMstep:  boolean determining whether an EM-step will be taken
#	  at each iteration
#	- ICMstep: boolean determining whether an ICM step will be taken
#	  at each iteration
#	- checkbnds: make sure that isotonization step does not wash out
#	  essential maximal antichains by using self-consistent bounds
#	- keepiter: boolean determining whether to keep the iteration
#	  states
#	- eps: maximal L1 distance between successive estimates before
#	  stopping iteration
#	- maxiter:  maximal number of iterations to perform before
#	  stopping
# Returned values:  object containing
#	- sigma:  NPMLE of the survival function on the maximal
#	  antichains
#	- weights:  diagonal of the likelihood function's second
#	  derivative
#	- lastchange:  vector of differences between the last two
#	  iterations
#	- numiter:  total number of iterations performed
#	- iter: (only present if keepiter is true) states of sigma during
#	  the iteration

VEMICMmac <- function(A,EMstep=TRUE,VEMstep=TRUE, ICMstep=TRUE, keepiter=FALSE, tol=1e-7,
                      tolbis=1e-7,maxiter=1000)
{
    if (!VEMstep && !ICMstep) {
        print("One of EMstep or ICMstep must be true.")
        return(NULL)
    }
    Meps<-.Machine$double.eps
    m<-dim(A)[1]
    n<-dim(A)[2]
    tA<-t(A)
    if (m==1) {
        ret<-NULL
        ret$sigma<-1
        ret$weights<-n
        ret$lastchange<-0
        ret$numiter<-0
        return(ret)
    }

    WW<-matrix(0,m,n)
    for (i in 1:(m-1))
        WW[i,]<-A[i,]-A[i+1,]
    WW[m,]<-A[m,]
    sigma<-cumsum(apply(A,1,sum)/sum(A))
    numiter<-0
    oldsigma<-rep(-1,m)
    if (keepiter) iter<-sigma
    while (max(abs(oldsigma-sigma))>tol && numiter<=maxiter) {
        oldsigma<-sigma
        if (VEMstep) {
            pvec<-diff(c(0,sigma))
            Linv<-1/(tA%*%pvec)
            dvec<-A%*%Linv
            mind<-min(dvec[pvec>0])
            minj<-match(mind,dvec)
            maxd<-max(dvec)
            maxj<-match(maxd,dvec)
            ndir<-rep(0,m)
            ndir[minj]<- -1
            ndir[maxj]<-1
            ndir<-ndir*pvec[minj]
            pold<-pvec
            pvec<-Bisect(tA,pvec,ndir,Meps,tolbis=tolbis)
            sigma<-cumsum(pvec)
            if (keepiter) iter<-rbind(iter,sigma)
        }
        if (ICMstep) {
            Wps<-1/(t(WW)%*%sigma)
            weights<-(abs(WW)%*%Wps^2)
            increment<-as.vector((WW%*%Wps)/weights)
            sigma<-sigma+increment
            sigma[m]<-1
            if (keepiter) iter<-rbind(iter,sigma)
            temp<-PMGA(sigma[-m],weights[-m])
            poolnum<-c(0,cumsum(temp$poolnum))
            for (i in 2:length(poolnum))
                for (j in (poolnum[i-1]+1):poolnum[i])
                    sigma[j]<-temp$est[i-1]
            if (keepiter) iter<-rbind(iter,sigma)
            # Implementing Jongbloed's correction through bisection
            temp<-c(0,sigma)
            pvec<-diff(c(0,oldsigma))
            ndir<-diff(c(0,temp[2:(m+1)]))-pvec
            pvec<-Bisect(tA,pvec,ndir,Meps,tolbis=tolbis)
            sigma<-cumsum(pvec)
            if (keepiter) iter<-rbind(iter,sigma)
        }
        numiter<-numiter+1
    }
    if (i==maxiter)
        warning("VEM/ICM may have failed to converge.")
    pf<-diff(c(0,sigma))
    ret<-list(sigma=sigma,pf=pf,lval=sum(log(t(A)%*%pf)),
              weights=as.vector(weights),lastchange=sigma-oldsigma,
              numiter=i,eps=tol)
    if (keepiter) {
        if (VEMstep && ICMstep)
            dimnames(iter)<-list(c("Seed",rep(c("VEM","Fisher",
                                                "PMGA","Bisect"),
                                              numiter)),NULL)
        else if (EMstep)
            dimnames(iter)<-list(rep("VEM",numiter+1),NULL)
        else
            dimnames(iter)<-list(c("Seed",rep(c("Fisher","PMGA"),
                                              numiter)),NULL)
        ret$iter<-iter
    }
    ret
}


#The EM algorithm for interval censored data

EM<-function(A, pvec, maxiter = 500, tol = 1e-12)
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





