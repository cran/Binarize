\name{plot,TrinarizationResult-method}
\alias{plot,TrinarizationResult-method}
\alias{plot,TrinarizationResult,ANY-method}
\alias{plot,numeric,TrinarizationResult-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of trinarization results.
}
\description{
	Visualizes a trinarization as a ray or a two-dimensional plot.
}
\usage{
\S4method{plot}{TrinarizationResult,ANY}(x,
	twoDimensional=FALSE,
	showLegend=TRUE,
	showThreshold=TRUE,
	...)
\S4method{plot}{numeric,TrinarizationResult}(x,
	y,
	showLegend=TRUE,
	showThreshold=TRUE,
	...)
}
\arguments{
	\item{x}{
		If \code{y} is supplied, this is a vector of x coordinates for the trinarization values in\code{y}, which are plotted on the y axis. If \code{y} is not supplied, this is object of class \code{\linkS4class{TrinarizationResult}} containing the trinarized values to visualize.
	}
	\item{y}{
		If \code{x} is a vector of x coordinates, this is object of class \code{\linkS4class{TrinarizationResult}} containing the trinarized values to visualize.
	}
	\item{twoDimensional}{
		Specifies whether the trinarization is depicted as a ray or as a two-dimensional curve (see details).
	}
	\item{showLegend}{
		If set to \code{true}, a legend is included in the plot.
	}
	\item{showThreshold}{
		If set to \code{true}, the trinarization thresholds are depicted as a horizontal or vertical lines (depending on \code{twoDimensional}).
	}
	\item{...}{
		Further graphical parameters to be passed to \code{\link{plot}}. The parameters \code{col} and \code{pch} can be supplied in different ways: If supplied as vectors of size 3, the first value corresponds to a 0 in the trinarization, the second value corresponds to a 1, and the third value corresponds to a 2. \code{col} can also have length 4, in which case the fourth entry is the color of the threshold line.
		If \code{col} or \code{pch} have the size of the input vector, the corresponding colors and symbols are assigned to the data points. 
	}
}
\details{
The function comprises two different plots: If \code{twoDimensional = TRUE}, the positions in the input vector are aligned with the x axis, and the y axis corresponds to the values. The trinarization thresholds are shown as a horizontal lines, and the trinarization is indicated by three different symbols.

If \code{twoDimensional = FALSE}, the trinarized values are aligned with a one-dimensional ray, and the separating thresholds are depicted as a vertical lines.
}

\seealso{
\code{\link{plot}}, 
\code{\linkS4class{TrinarizationResult}}
}
\examples{
# plot a binarization in one and two dimensions
res <- TASC(iris[,"Petal.Length"])
plot(res)
plot(res, twoDimensional = TRUE)
plot(res, twoDimensional = TRUE, 
	 pch = c("x", "+"), 
	 col = c("red", "black", "royalblue", "green"), 
	 lty = 4, lwd = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
