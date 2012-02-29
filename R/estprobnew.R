## file: estprobnew.R
## function: estprobnew
## Usage: estimaing the probability that a randomly sampled n-th respondent
## states a new WTP-interval which differs from any observed in the sample
## of size n-1. The estimated probability is decreasing function of n.
## The estimator is obtained by the conditional averaging 
## Date created: 2010-6-15
## input: data -- cases of self-selected intervals
##        size -- random permutation size
## output: object of class 'iwtp' 
##

##th = c(2, 3, 4, 1, 1, 1, 1, 5, 5, 1, 6, 2, 3, 5, 39, 1, 1, 1, 11, 1, 
##   3, 1, 1, 1, 2, 3, 69, 3, 3, 1, 1, 8, 1, 23, 3, 4, 5, 1, 2, 3, 1, 
##   1, 1, 1, 4, 1)

estprobnew <- function(th, size=500) {
	n <- length(th)
	dh <- c()
	for(i in 1:n) {
		 dh <- c(dh, rep(i, th[i]))
	}
	#dh
	nd <- length(dh)

	ch <- 0
	set.seed(5743)
	for(s in 1:size) {
		rp <- sample(1:nd, replace=FALSE)
		#rp

		dp <- dh[rp]
		#dp

		ld <- dp[nd]
		#ld

		dd <- dp[1:(nd-1)]
		#dd

		isnew <- ifelse(length(which(dd==ld))>0,0,1)

		ch <- ch + isnew
	}

	return(ch/size)
}