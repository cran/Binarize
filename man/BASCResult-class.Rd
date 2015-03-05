\name{BASCResult-class}
\Rdversion{1.1}
\docType{class}
\alias{BASCResult-class}
\alias{plotStepFunctions,BASCResult-method}
\alias{print,BASCResult-method}
\alias{show,BASCResult-method}

\title{Class "BASCResult"}
\description{
    A specialized class storing the results of a call to \code{\link{binarize.BASC}}.
}
\section{Objects of this class}{
    Objects of this class shouldn't be created directly. They are created implicitly by a call to \code{\link{binarize.BASC}}.
}
\section{Slots}{
  \describe{
    \item{\code{p.value}:}{The p-value of the statistical test for reliability of the binarization.}
    \item{\code{intermediateSteps}:}{A matrix specifying the optimal step functions from which the binarization was calculated. The number of rows corresponds to the number of step functions, and the number of columns is determined by the length of the input vector minus 2 (that is, the length of the step function corresponding to the input vector). From the first to the last row, the number of steps increases. The non-zero entries of the matrix represent the locations of the steps. Step functions with fewer steps than the input step function have entries set to zero.}
    \item{\code{intermediateHeights}:}{A matrix giving the jump heights of the steps supplied in \code{intermediateSteps}.}
    \item{\code{intermediateStrongestSteps}:}{A vector with one entry for each step function (row) in \code{intermediateSteps}. The entries specify the location of the strongest step for each of the functions.}
   \item{\code{originalMeasurements}:}{A numeric vector storing the input measurements.}
    \item{\code{binarizedMeasurements}:}{An integer vector of binarized values (0 or 1) corresponding to the original measurements.}
    \item{\code{threshold}:}{The threshold that separates 0 and 1.}
    \item{\code{method}:}{A string describing the binarization method that yielded the result.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{BinarizationResult}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\link{plotStepFunctions}}{\code{signature(x = "BASCResult")}: Plot the intermediate optimal step functions used to determine the threshold.}
    \item{print}{\code{signature(x = "BASCResult")}: Print a summary of the binarization.}
    \item{show}{\code{signature(object = "BASCResult")}: ... }
	 }
}

\seealso{
\code{\link{binarize.BASC}}, 
\code{\linkS4class{BinarizationResult}}
}
\keyword{classes}
