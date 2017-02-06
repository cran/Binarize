\name{binarizeMatrix}
\alias{binarizeMatrix}
\title{
Utility function to binarize a matrix of measurements
}
\description{
Binarizes a matrix of measurements all at once, and returns the binarized vectors as well as the binarization thresholds and the p-values. 
}
\usage{
binarizeMatrix(mat, 
               method = c("BASCA", "BASCB", "kMeans"), 
               adjustment = "none", 
               ...)
}
\arguments{
  \item{mat}{
A n x m matrix comprising m raw measurements of n features.
}
  \item{method}{
The binarization algorithm to be used. \code{method="BASCA"} calls \code{\link{binarize.BASC}} with \code{method="A"}. \code{method="BASCB"} calls \code{\link{binarize.BASC}} with \code{method="B"}. \code{method="kMeans"} calls \code{\link{binarize.kMeans}}.
}
  \item{adjustment}{
Specifies an optional adjustment for multiple testing that is applied to the p-values (see \code{\link{p.adjust}} for possible values). By default, no adjustment is applied.
}
  \item{\dots}{
Further parameters that are passed to the respective binarization methods (\code{\link{binarize.BASC}}  or \code{method="kMeans"}).
}
}
\value{
A n x (m+2) matrix of binarized measurements. Here, the first m columns correspond to the binarized measurements. The m+1-st column comprises the binarization thresholds for the features, and the m+2-nd column contains the p-values.
}

\seealso{
\code{\link{binarize.BASC}}, \code{\link{binarize.kMeans}}, \code{\link{p.adjust}}
}
\examples{
bin <- binarizeMatrix(t(iris[,1:4]))
print(bin)
}
\keyword{ binarization }
\keyword{ binarize }
\keyword{ BASC }
\keyword{ k-means }
\keyword{ matrix }
\keyword{ multiple scales }
