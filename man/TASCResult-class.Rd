\name{TASCResult-class}
\Rdversion{1.1}
\docType{class}
\alias{TASCResult-class}
\alias{plotStepFunctions,TASCResult-method}
\alias{print,TASCResult-method}
\alias{show,TASCResult-method}

\title{Class "TASCResult"}
\description{
	A specialized class storing the results of a call to \code{\link{TASC}}.
}
\section{Objects of this class}{
	Objects of this class shouldn't be created directly. They are created implicitly by a call to \code{\link{TASC}}.
}
\section{Slots}{
	\describe{
		\item{\code{p.value}:}{The p-value of the statistical test for reliability of the trinarization.}
		\item{\code{intermediateSteps}:}{A matrix specifying the optimal step functions from which the trinarization was calculated. The number of rows corresponds to the number of step functions, and the number of columns is determined by the length of the input vector minus 2 (that is, the length of the step function corresponding to the input vector). From the first to the last row, the number of steps increases. The non-zero entries of the matrix represent the locations of the steps. Step functions with fewer steps than the input step function have entries set to zero.}
		\item{\code{intermediateHeights1}:}{A matrix giving the jump heights of the steps supplied in \code{intermediateSteps} for the first threshold.}
		\item{\code{intermediateHeights2}:}{A matrix giving the jump heights of the steps supplied in \code{intermediateSteps} for the second threshold.}
		\item{\code{intermediateStrongestSteps}:}{A matrix with one row for each step function (row) in \code{intermediateSteps}. The entries specify the location of the two strongest steps for each of the functions.}
		\item{\code{originalMeasurements}:}{A numeric vector storing the input measurements.}
		\item{\code{trinarizedMeasurements}:}{An integer vector of trinarized values (0, 1 or 2) corresponding to the original measurements.}
		\item{\code{threshold1}:}{The threshold that separates 0 from 1.}
		\item{\code{threshold2}:}{The threshold that separates 1 from 2.}
		\item{\code{method}:}{A string describing the trinarization method that yielded the result.}
	}
}
\section{Extends}{
Class \code{"\linkS4class{TrinarizationResult}"}, directly.
}
\section{Methods}{
	\describe{
		\item{\link{plotStepFunctions}}{\code{signature(x = "TASCResult")}: Plot the intermediate optimal step functions used to determine the thresholds.}
		\item{print}{\code{signature(x = "TASCResult")}: Print a summary of the trinarization.}
		\item{show}{\code{signature(object = "TASCResult")}: ... }
	}
}

\seealso{
\code{\link{TASC}}, 
\code{\linkS4class{TrinarizationResult}}
}
\keyword{classes}
