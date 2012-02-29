\encoding{CP1252}
\name{resam}
\Rdversion{1.1}
\alias{resam}
\alias{QQplot}
\alias{QQplot.iwtp.resam}
\alias{plot.iwtp.resam.wb}
\alias{print.iwtp.resam.wb}
\alias{print.iwtp.resam.wemix}
\title{
Resampling methods for estimation accuracy of statistical inferecnes
}
\description{
Function \code{resam} carries out re-samplings for a given number times. For
each copy of the re-sampled data, the survival function is estimated. 
Deviations of ML-esimators from mean WTP are calculated. The obtained
deviations can be used to draw normal Quantile-Quantile (Q-Q) plot.
}
\usage{
resam(obj,size=10,seed=1234,plot=FALSE,
			bounds=list(),limits=list())
\method{QQplot}{iwtp.resam}(x,x2=NULL,xlim=NULL,ylim=NULL,
			xlab="",ylab="",col="red",lty=1,lwd=1,asp=0.8,\dots)
\method{plot}{iwtp.resam.wb}(x,xlab="",ylab="",\dots)
\method{print}{iwtp.resam.wb}(x,\dots)
\method{print}{iwtp.resam.wemix}(x,\dots)
}

\arguments{
  \item{obj}{
object of class 'iwtp'
}
  \item{size}{
integer, sampling size
}
  \item{seed}{
seed for random number
}
  \item{plot}{
logical, whether to present a scatter plot for the estimatetd
  parameters: a, b for the case of Weibull distribution
}
  \item{bounds}{
a list, specify the lower and upper bounds and the initial values for the parameters, See function \code{iwtp}
}

  \item{limits}{
a list,  specify the lower and upper bounds and the initial values for the parameters of the function for modeling behaviour, See function \code{iwtp}
}
\item{x,x2}{
output from function \code{resam}
  }
  \item{xlim,ylim}{
     numerical vectors of length 2, giving the x and y coordinates ranges
  }
  \item{xlab}{
   character, title of x-axis
  }
  \item{ylab}{
   character, title of y-axis
  }
  \item{col}{
  line color, it can be a vector when more than one line are
  drawn, See function \code{iwtp}
  }
  \item{lty}{
  line type, it can be a vector when more than one line are drawn, See
  function \code{iwtp}
  }
  \item{lwd}{
  line width, it can be a vector when more than one line are drawn, See
  function \code{iwtp}
  }
  \item{asp}{
    numerical, giving the aspect ratio y/x, See function \code{asp}
  }
  \item{\dots}{
further graphical parameters as in \code{par}
  }
}
\details{
Argument \code{obj} must be a object of class 'iwtp', it is an output of function \code{iwpt}.
Argument \code{size} gives the number of resamplings. Freq.redo tells the function whether to re-estimate
the probability associated with each of the division intervals.
}
\value{
Function \code{resam} returns a object of class 'iwtp.resam.wb' or 'iwtp.resam.wemix', depending on the obj. 
The object is a LIST, which consists of:
 \item{par}{a matrix with the estimated parameters for each copy of resamplings.}
 \item{sf}{survival function, a copy from argument obj}
 \item{size}{sample size, a copy of argument size}
 \item{dev.mwtp}{a numerical vector, the deviations from estimated mean WTP}
 \item{run.time}{the run time taken by the resampling}
}
\references{
Belyaev, Yuri and Bengt Kriström. 2010. Approach to Analysis of Self-Selected Interval
Data. Technical report SLU, Department of Forest Economics 90183 Umeå, Sweden.
}
\author{
Wenchao Zhou
}

\seealso{
\code{\link{iwtp}}
}
\examples{

data(IL)

## the Weibull survival function
wtp.wb <- iwtp(IL,"weibull")

## Resample without re-estimation of frequencies
res.wb <- resam(wtp.wb, size=10,plot=FALSE)

## scatter plot
plot(res.wb)

## show QQ plot
QQplot(res.wb)

## the mixed Weibull/exponential survival function
wtp.wemix <- iwtp(IL, dist="wemix")
##
res.wemix <- resam(wtp.wemix, size=10, plot=FALSE)
QQplot(res.wemix)
QQplot(res.wb,res.wemix,col=c("red","blue"))
}

\keyword{survival}

