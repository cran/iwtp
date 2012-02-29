\encoding{CP1252}
\name{iwtp}
\Rdversion{1.1}
\alias{iwtp}
\alias{plot.iwtp}
\alias{print.iwtp}

\title{
Estimation of Survival Function (Weibull or mixed Weibull)
}

\description{
Estimating the Weibulll or the mixed Weibull and exponential distribution function
}

\usage{
iwtp(data,dist="weibull",bm.type=5,plot=FALSE,bounds=list(),limits=list())
\S3method{plot}{iwtp}(x,x2=NULL,xlim=NULL,xlab="",col="red",lty=1,lwd=1,asp=0.6,\dots)
\S3method{print}{iwtp}(x,\dots)
}

\arguments{
  \item{data}{
a data.frame
}
	\item{bm.type}{
	numeric number, indicating the type of behavior of respondents selecting intervals. It can be 1, 2, 3, 4 and 5
	}
  \item{dist}{
character, type of survival function. It can be "weibull" or "wemix", See Details
}
  \item{plot}{
logical, whether to plot the empirical survival function
}
  \item{bounds}{
a list, specify the lower and upper bounds and the initial values for the parameters, which will be used by the optimizer, See \code{optim}
}

  \item{limits}{
a list, specify the lower and upper bounds and the initial values for the function of modelling behavior, which will be
  used by the optimizer, See \code{optim}
}
\item{x,x2}{
the output from function \code{iwtp}
  }
  \item{xlim}{
    numeric vector with length 2, giving the x-axis range
  }
  \item{xlab}{
    character, title of x-axis
  }
  \item{col}{
  line color, it can be a vector when more than one line are
  drawn, See \code{par}
  }
  \item{lty}{
  line type, it can be a vector when more than one line are drawn, See \code{par}
  }
  \item{lwd}{
  line width, it can be a vector when more than one line are drawn, See \code{par}
  }
  \item{asp}{
	numerical, giving the aspect ratio y/x, See \code{par}.
  }
  \item{\dots}{
further graphical parameters as in \code{par}.
  }
}
\details{
Argument \code{data} must be a data.frame. \code{data} should contain at least 4 columns. 
Its first column is the case number, column 2 is the left-end value of stated intervals, 
column 2 is the right-end value of stated intervals, and column 4 is the times of 
the interval being stated.

Argument \code{dist} can be "weibull" or "wemix". If "weibull" is chosen, a
Weibull survival function will be estimated; If "wemix" is chosen, a
mixed Weilbull/expoential survival function will be estimated.

}
\value{
Function \code{iwtp} returns an object of class 'iwtp'. An iwtp object is a LIST with the
following three components:
  \item{sf}{parameter estimates for the assigned survival function}
  \item{esf}{a data list of the empirical survival function based on the right-end value of stated intervals}
  \item{intv}{an object of class 'iwtp.interv', containing information about division intervals}
}
\references{
Belyaev, Yuri and Bengt Kriström. 2010. Approach to Analysis of Self-Selected Interval
Data. Technical report SLU, Department of Forest Economics 90183 Umeå, Sweden.
}
\author{
Wenchao Zhou
}


\seealso{
\code{\link{optim}}
}
\examples{

data(IL)

# survival function: the Weibull distribution
wtp.wb <- iwtp(IL,dist="weibull",plot=TRUE)
wtp.wb

#surival function: the mixed Weibull and exponential distribution
wtp.wemix <- iwtp(IL, dist="wemix", plot=TRUE)
wtp.wemix


# setting initial values for parameters
bounds <- list(
		a=c(lower=70,upper=300,init=75),
		b=c(lower=1,upper=4.5,init=1.4))
wtp.wb <- iwtp(IL, "weibull",plot=FALSE,bounds=bounds)

# plot empirical survival function
esfPlot(wtp.wb,col="red",lwd=2)

# plot empirical s.f and estimated s.fs in a single figure
plot(wtp.wb,wtp.wemix,col=c("red","green"),lty=c(1,2))
}

\keyword{survival}

