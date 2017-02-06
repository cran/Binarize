\name{plot,BinarizationResult-method}
\alias{plot,BinarizationResult-method}
\alias{plot,BinarizationResult,ANY-method}
\alias{plot,numeric,BinarizationResult-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of binarization results.
}
\description{
  Visualizes a binarization as a ray or a two-dimensional plot. 
}
\usage{
\S4method{plot}{BinarizationResult,ANY}(x,
     twoDimensional=FALSE, 
     showLegend=TRUE, 
     showThreshold=TRUE, 
     ...)
\S4method{plot}{numeric,BinarizationResult}(x, 
                                            y, 
                                            showLegend=TRUE, 
                                            showThreshold=TRUE, 
                                            ...)     
}
\arguments{
    \item{x}{
 If \code{y} is supplied, this is a vector of x coordinates for the binarization values in\code{y}, which are plotted on the y axis. If \code{y} is not supplied, this is object of class \code{\linkS4class{BinarizationResult}} containing the binarized values to visualize.
    }
    
    \item{y}{
If \code{x} is a vector of x coordinates, this is object of class \code{\linkS4class{BinarizationResult}} containing the binarized values to visualize.    
    }    
    \item{twoDimensional}{
 Specifies whether the binarization is depicted as a ray or as a two-dimensional curve (see details).
    }
    \item{showLegend}{
 If set to \code{true}, a legend is included in the plot.
    }
    \item{showThreshold}{
 If set to \code{true}, the binarization threshold is depicted as a horizontal or vertical line (depending on \code{twoDimensional}).
    }
    \item{...}{
 Further graphical parameters to be passed to \code{\link{plot}}. The parameters \code{col} and \code{pch} can be supplied in different ways: If supplied as vectors of size 2, the first value corresponds to a 0 in the binarization, and the second value corresponds to a 1 in the binarization. \code{col} can also have length 3, in which case the third entry is the color of the threshold line. 
If \code{col} or \code{pch} have the size of the input vector, the corresponding colors and symbols are assigned to the data points. 
}
}
\details{
The function comprises two different plots: If \code{twoDimensional = TRUE}, the positions in the input vector are aligned with the x axis, and the y axis corresponds to the values. The binarization threshold is shown as a horizontal line, and the binarization is indicated by two different symbols. 
  
If \code{twoDimensional = FALSE}, the binarized values are aligned with a one-dimensional ray, and the separating threshold is depicted as a vertical line. 
}

\seealso{
\code{\link{plot}}, 
\code{\linkS4class{BinarizationResult}}
}
\examples{
# plot a binarization in one and two dimensions
res <- binarize.BASC(iris[,"Petal.Length"], method="A")
plot(res)
plot(res, twoDimensional = TRUE)
plot(res, twoDimensional = TRUE, 
     pch = c("x", "+"), 
     col = c("red", "black", "royalblue"), 
     lty = 4, lwd = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
