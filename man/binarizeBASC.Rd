\name{binarize.BASC}
\alias{binarize.BASC}
\title{
Binarization Across Multiple Scales
}
\description{Binarizes real-valued data using the multiscale BASC methods.}
\usage{
binarize.BASC(vect, 
              method = c("A","B"), 
              tau = 0.01, 
              numberOfSamples = 999, 
              sigma = seq(0.1, 20, by=.1),
              na.rm=FALSE)
}

\arguments{
  \item{method}{
    Chooses the BASC method to use (see details), i.e. either  "A" or "B". 
  }
  \item{vect}{
    A real-valued vector of data to binarize.
  }
  \item{tau}{
    This parameter adjusts the sensitivity and the specificity of the statistical testing procedure that rates the quality of the binarization. Defaults to 0.01.
  }
  \item{numberOfSamples}{
    The number of samples for the bootstrap test. Defaults to 999.
  }
  \item{sigma}{
    If \code{method="B"}, this specifies a vector of different sigma values for the convolutions with the Bessel function. Ignored for \code{method="A"}.
  }
  
  \item{na.rm}{
  If set to \code{TRUE}, \code{NA} values are removed from the input. Otherwise, binarization will fail in the presence of \code{NA} values.
  }
}
\details{
The two BASC methods can be subdivided into three steps:
\describe{	
	\item{Compute a series of step functions:}{An initial step function is obtained by rearranging the original time series measurements in increasing order. Then, step functions with fewer discontinuities are calculated. BASC A calculates these step functions in such a way that each minimizes the Euclidean distance to the initial step function. BASC B obtains step functions from smoothened versions of the input function in a scale-space manner.}
	\item{Find strongest discontinuity in each step function:}{A strong discontinuity is a high jump size (derivative) in combination with a low approximation error.}
	\item{Estimate location and variation of the strongest discontinuities:}{Based on these estimates, data values can be excluded from further analyses.
}
}
}
\value{
    Returns an object of class \code{\linkS4class{BASCResult}}.
}
\references{
M. Hopfensitz, C. Müssel, C. Wawra, M. Maucher, M. Kuehl, H. Neumann, and H. A. Kestler. Multiscale Binarization of Gene Expression Data for Reconstructing Boolean Networks. IEEE/ACM Transactions on Computational Biology and Bioinformatics 9(2):487-498, 2012.).
}

\seealso{
    \code{\linkS4class{BinarizationResult}}, 
    \code{\linkS4class{BASCResult}}
}
\examples{
par(mfrow=c(2,1))
result <- binarize.BASC(iris[,"Petal.Length"], method="A", tau=0.15)
print(result)
plot(result)

result <- binarize.BASC(iris[,"Petal.Length"], method="B", tau=0.15)
print(result)
plot(result)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ binarization }
\keyword{ binarize }
\keyword{ BASC }
\keyword{ multiple scales }
