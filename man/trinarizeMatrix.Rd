\name{trinarizeMatrix}
\alias{trinarizeMatrix}
\title{
Utility function to trinarize a matrix of measurements
}
\description{
Trinarizes a matrix of measurements all at once, and returns the trinarized vectors as well as the trinarization thresholds and the p-values. 
}
\usage{
trinarizeMatrix(mat, 
               method = c("TASCA", "TASCB"), 
               adjustment = "none", 
               ...)
}
\arguments{
\item{mat}{
	A n x m matrix comprising m raw measurements of n features.
}
\item{method}{
	The trinarization algorithm to be used. \code{method="TASCA"} calls \code{\link{TASC}} with \code{method="A"}. \code{method="TASCB"} calls \code{\link{TASC}} with \code{method="B"}.
}
\item{adjustment}{
	Specifies an optional adjustment for multiple testing that is applied to the p-values (see \code{\link{p.adjust}} for possible values). By default, no adjustment is applied.
}
\item{\dots}{
	Further parameters that are passed to the respective trinarization methods (\code{\link{TASC}}).
}
}
\value{
	A n x (m+3) matrix of trinarized measurements. Here, the first m columns correspond to the trinarized measurements. The m+1-st and the m+2-st column comprises the trinarization thresholds for the features, and the m+3-nd column contains the p-values.
}

\seealso{
\code{\link{TASC}}, \code{\link{p.adjust}}
}
\examples{
tri <- trinarizeMatrix(t(iris[,1:4]))
print(tri)
}
\keyword{ trinarization }
\keyword{ trinarize }
\keyword{ TASC }
\keyword{ matrix }
\keyword{ multiple scales }
