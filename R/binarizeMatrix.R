binarizeMatrix <- function(mat, method=c("BASCA","BASCB","kMeans"), adjustment="none", ...)
{
  binFunc <- switch(match.arg(method, c("BASCA","BASCB","kMeans")),
    "BASCA" = function(x, ...)
    {
      bin <- binarize.BASC(x, method="A", ...)
      return(c(bin@binarizedMeasurements, 
               list(bin@threshold, bin@p.value)))
    },
    "BASCB" = function(x, ...)
    {
      bin <- binarize.BASC(x, method="B", ...)
      return(c(bin@binarizedMeasurements, 
              list(bin@threshold, bin@p.value)))
    },
    "kMeans" = function(x, ...)
    {
      bin <- binarize.kMeans(x, ...)
      return(c(bin@binarizedMeasurements, 
               list(bin@threshold, bin@p.value)))
    }
  )

  bin <- do.call("rbind.data.frame", apply(mat, 1, binFunc, ...))
  
  if (!is.null(colnames(mat)))
    colnames(bin) <- c(colnames(mat), "threshold", "p.value")
  else
    colnames(bin) <- c(paste("V",seq_len(ncol(mat)),sep=""), 
                       "threshold", "p.value")

                         
  bin[,"p.value"] <- p.adjust(bin[,"p.value"], method=adjustment)
  return(bin)
}
