# TODO: Add comment
# 
# Author: zhouw
###############################################################################
#turnbull style interval data weibull

turnbull.weibull <- function(data)
{
	ok <- require(survival)
	if(!ok)
	{
		stop("Package survival must be installed!", call. = FALSE)
	}
	left <- data[,2]
	# the left value cann't be 0
	ndx <- which(left==0)
	left[ndx] <- 0.0001+left[ndx]
	# call function Surv from package survival
	s.me=Surv(left, data[,3], type="interval2" )
	res.me=survreg(s.me~1)
	#coef(res.me)
	scale.me=res.me$scale
	mean.me=gamma(1+scale.me)*exp(coef(res.me))
	mean.me
}