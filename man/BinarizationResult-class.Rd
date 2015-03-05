\name{BinarizationResult-class}
\Rdversion{1.1}
\docType{class}
\alias{BinarizationResult-class}
\alias{plotBinarization,BinarizationResult-method}
\alias{print,BinarizationResult-method}
\alias{show,BinarizationResult-method}

\title{Class "BinarizationResult"}
\description{
    This is the base class for objects that store the results of a binarization algorithm. It defines the slots and methods that the results of all algorithms share.
}
\section{Objects of this class}{
    Objects of this class shouldn't be created directly. They are created implicitly by a call to one of the binarizeation algorithms.
}
\section{Slots}{
  \describe{
    \item{\code{originalMeasurements}:}{A numeric vector storing the input measurements.}
    \item{\code{binarizedMeasurements}:}{An integer vector of binarized values (0 or 1) corresponding to the original measurements.}
    \item{\code{threshold}:}{The threshold that separates 0 and 1.}
    \item{\code{method}:}{A string describing the binarization method that yielded the result.}
    \item{\code{p.value}:}{The p-value obtained by a test for validity of the binarization (e.g. BASC bootstrap test, Hartigan's dip test for k-means binarization, scan statistic p-value for best window. If no test was performed, this is \code{NA}.}
  }
}
\section{Methods}{
  \describe{
    \item{\link[=plot,BinarizationResult-method]{plot}}{\code{signature(x = "BinarizationResult")}: Plot the binarization and the threshold.}
    \item{print}{\code{signature(x = "BinarizationResult")}: Print a summary of the binarization.}
    \item{show}{\code{signature(object = "BinarizationResult")}: ... }
	 }
}

\seealso{
\code{\link{binarize.BASC}}, 
%\code{\link{binarize.scanStatistic}}, 
%\code{\link{binarize.edgeDetector}}, 
\code{\link{binarize.kMeans}}, 
\code{\linkS4class{BASCResult}}, 
%\code{\linkS4class{ScanStatisticResult}}
}
\keyword{classes}
