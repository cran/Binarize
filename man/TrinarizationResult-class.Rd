\name{TrinarizationResult-class}
\Rdversion{1.1}
\docType{class}
\alias{TrinarizationResult-class}
\alias{plotTrinarization,TrinarizationResult-method}
\alias{print,TrinarizationResult-method}
\alias{show,TrinarizationResult-method}

\title{Class "TrinarizationResult"}
\description{
    This is the base class for objects that store the results of a trinarization algorithm. It defines the slots and methods that the results of all algorithms share.
}
\section{Objects of this class}{
    Objects of this class shouldn't be created directly. They are created implicitly by a call to one of the trinarizeation algorithms.
}
\section{Slots}{
  \describe{
    \item{\code{originalMeasurements}:}{A numeric vector storing the input measurements.}
    \item{\code{trinarizedMeasurements}:}{An integer vector of trinarized values (0 or 1 or 2) corresponding to the original measurements.}
    \item{\code{threshold1}:}{The threshold that separates 0 and 1.}
    \item{\code{threshold2}:}{The threshold that separates 1 and 2.}
    \item{\code{method}:}{A string describing the trinarization method that yielded the result.}
    \item{\code{p.value}:}{The p-value obtained by a test for validity of the trinarization (e.g. TASC bootstrap test). If no test was performed, this is \code{NA}.}
  }
}
\section{Methods}{
  \describe{
    \item{\link[=plot,BinarizationResult-method]{plot}}{\code{signature(x = "TrinarizationResult")}: Plot the trinarization and the thresholds.}
    \item{print}{\code{signature(x = "TrinarizationResult")}: Print a summary of the trinarization.}
    \item{show}{\code{signature(object = "TrinarizationResult")}: ... }
	 }
}

\seealso{
\code{\link{TASC}}, 
%\code{\link{binarize.scanStatistic}}, 
%\code{\link{binarize.edgeDetector}}, 
%\code{\link{binarize.kMeans}}, 
\code{\linkS4class{TASCResult}}, 
%\code{\linkS4class{ScanStatisticResult}}
}
\keyword{classes}
