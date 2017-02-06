\name{plotStepFunctions}
\alias{plotStepFunctions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot all step functions for BASC or TASC
}
\description{
A specialized visualization that plots all the optimal step functions computed by the BASC algorithms or TASC.
}
\usage{
plotStepFunctions(x, 
                  showLegend=TRUE, 
                  connected=FALSE, 
                  withOriginal=TRUE, 
                  ...)
}
\arguments{
  \item{x}{
   A binarization (or trinarisation) result object of class \code{\linkS4class{BASCResult}} (\code{\linkS4class{TASCResult}}).
  }
  \item{showLegend}{
   If \code{TRUE}, a legend is included in the plot.
  }
  \item{connected}{
   If \code{TRUE}, the single steps of the step functions are connected by lines.
  }
  \item{withOriginal}{
   If set to \code{TRUE}, the original step function (i.e. the sorted input vector) is included in the plot.
  }
  \item{...}{
   Additional graphical parameters to be passed to \code{\link{plot}}.
  }
}

\seealso{
\code{\linkS4class{BASCResult}},
\code{\link{binarize.BASC}},
\code{\linkS4class{TASCResult}},
\code{\link{TASC}} 
}
\examples{
result <- binarize.BASC(iris[,"Petal.Width"], 
                        method="B")
plotStepFunctions(result)

result <- TASC(iris[,"Petal.Width"])
plotStepFunctions(result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
