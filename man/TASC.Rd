\name{TASC}
\alias{TASC}
\title{
Trinarization Across Multiple Scales
}
\description{Trinarizes real-valued data using the multiscale TASC method.}
\usage{
TASC(vect, 
	method = c("A","B"), 
	tau = 0.01, 
	numberOfSamples = 999, 
	sigma = seq(0.1, 20, by=.1),
	na.rm=FALSE,
	error = c("mean", "min"))
}

\arguments{
	\item{method}{
		Chooses the TASC method to use (see details), i.e. either  "A" or "B". 
	}
	\item{vect}{
	A real-valued vector of data to trinarize.
	}
	\item{tau}{
		This parameter adjusts the sensitivity and the specificity of the statistical testing procedure that rates the quality of the trinarization. Defaults to 0.01.
	}
	\item{numberOfSamples}{
		The number of samples for the bootstrap test. Defaults to 999.
	}
	\item{sigma}{
		If \code{method="B"}, this specifies a vector of different sigma values for the convolutions with the Bessel function. Ignored for \code{method="A"}.
	}
	\item{na.rm}{
		If set to \code{TRUE}, \code{NA} values are removed from the input. Otherwise, trinarization will fail in the presence of \code{NA} values.
	}
	\item{error}{
		Determines which error should be used for the data points between two thresholds, the "mean" error (default) to the thresholds or the "min" error.
	}
}
\details{
	The two TASC methods can be subdivided into three steps:
	\describe{	
		\item{Compute a series of step functions:}{An initial step function is obtained by rearranging the original time series measurements in increasing order. Then, step functions with fewer discontinuities are calculated. TASC A calculates these step functions in such a way that each minimizes the Euclidean distance to the initial step function. TASC B obtains step functions from smoothened versions of the input function in a scale-space manner.}
		\item{Find strongest discontinuities in each step function:}{A strong discontinuity is a high jump size (derivative) in combination with a low approximation error. For TASC a pair of strongest discontinuities is determined.}
		\item{Estimate location and variation of the strongest discontinuities:}{Based on these estimates, data values can be excluded from further analyses.
		}
	}
}
\value{
	Returns an object of class \code{\linkS4class{TASCResult}}.
}

\seealso{
	\code{\linkS4class{TrinarizationResult}}, 
	\code{\linkS4class{TASCResult}}
}
\examples{
par(mfrow=c(2,1))
result <- TASC(iris[,"Petal.Width"], method="A", tau=0.15)
print(result)
plot(result)

result <- TASC(iris[,"Petal.Width"], method="B", tau=0.15)
print(result)
plot(result)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ binarization }
\keyword{ binarize }
\keyword{ BASC }
\keyword{ multiple scales }
