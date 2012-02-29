\encoding{CP1252}
\name{iwtp.esf}
\Rdversion{1.1}
\alias{esfPlot}
\alias{iwtp.esf}
\alias{esfPlot.iwtp}
\title{
Empirical survival function
}
\description{
This function builds the empirical survival function for a given data list
}
\usage{
iwtp.esf(data)
\method{esfPlot}{iwtp}(x,x2=NULL,xlim=NULL,xlab="",col="red",lty=1,lwd=1,asp=0.6,\dots)
}

\arguments{
  \item{data}{
	data list, it must be a matrix with at least two columns. Coumn 1 is the point estimates, and comlumn 2
	is the times of that point observation being stated.
}
\item{x,x2}{
output from function \code{iwtp}
  }
  \item{xlim}{
    numerical vectors of length 2, giving the x ranges 
  }
  \item{xlab}{
  character,  title for x-axis
  }
  \item{col}{
  line color, it can be a vector when more than one line are drawn, See \code{par}
  }
  \item{lty}{
  line type, it can be a vector when more than one line are drawn, See \code{par}
  }
  \item{lwd}{
  line width, it can be a vector when more than one line are drawn, See \code{par}
  }
  \item{asp}{
	numerical, giving the aspect ratio y/x, See \code{par}
  }
  \item{\dots}{
further graphical parameters as in \code{par}. 
  }


}
\details{
	Argument data must be a matrix with at least two columns. The nonparametric mean will be calculated.
}
\value{
Function \code{iwtp.esf} returns a object of class 'iwtp.esf', which is a LIST.
\item{esf}{data list of the empirical survival function}
\item{mean}{the mean value of the empirical survival function}
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
esf <- iwtp.esf(IL)
# show mean wtp
esf$mean

# plot the empirical survival function
esfPlot(esf, col="blue",lty=1,lwd=2)
}

\keyword{ survival }

