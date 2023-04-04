\name{binarize.kMeans}
\alias{binarize.kMeans}
\title{
k-means Binarization
}
\description{
Binarizes a vector of real-valued data using the k-means clustering algorithm. The data is first split into 2 clusters.The values belonging to the cluster with the smaller centroid are set to 0, and the values belonging to the greater centroid are set to 1.
}
\usage{
binarize.kMeans(vect, 
                nstart=1, 
                iter.max=10,
                dip.test=TRUE,
                na.rm=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{vect}{
  A real-valued vector to be binarized (at least 3 measurements).
    }
    \item{nstart}{
  The number of restarts for k-means. See \code{\link{kmeans}} for details.}
    \item{iter.max}{
  The maximum number of iterations for k-means. See \code{\link{kmeans}} for details.}
    \item{dip.test}{If set to \code{TRUE}, Hartigan's dip test for unimodality is performed on \code{vect}, and its p-value is returned in the \code{pvalue} slot of the result. An insignificant test indicates that the data may not be binarizeable.}
    
      \item{na.rm}{
  If set to \code{TRUE}, \code{NA} values are removed from the input. Otherwise, binarization will fail in the presence of \code{NA} values.
  }
}
\value{
Returns an object of class \code{\linkS4class{BinarizationResult}}. 
}

\examples{
result <- binarize.kMeans(iris[,"Petal.Length"])

print(result)
plot(result, twoDimensional=TRUE)
}
\seealso{
\code{\link{kmeans}}, 
\code{\linkS4class{BinarizationResult}}, 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
