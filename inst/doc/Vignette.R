### R code from vignette source 'Vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Vignette.Rnw:38-39
###################################################
set.seed(13579)


###################################################
### code chunk number 2: Vignette.Rnw:156-157 (eval = FALSE)
###################################################
## install.packages("Binarize")


###################################################
### code chunk number 3: Vignette.Rnw:160-161
###################################################
library("Binarize")


###################################################
### code chunk number 4: Vignette.Rnw:174-175
###################################################
data(binarizationExample)


###################################################
### code chunk number 5: Vignette.Rnw:181-188
###################################################
pdf("density.pdf")
par(mar=c(2,2,1,1))
#plot(density(binarizationExample[1,]),main="")
#abline(v=mean(binarizationExample[1,]), lty="dashed")
plot(function(x)dnorm(x,mean=0,sd=1)+dnorm(x,mean=10,sd=1),xlim=c(-5,15),main="")
abline(v=5, lty="dashed")
dev.off()


###################################################
### code chunk number 6: Vignette.Rnw:202-204
###################################################
bin <- binarize.kMeans(binarizationExample[1,])
print(bin)


###################################################
### code chunk number 7: Vignette.Rnw:211-212
###################################################
print(bin@binarizedMeasurements)


###################################################
### code chunk number 8: Vignette.Rnw:216-217 (eval = FALSE)
###################################################
## plot(bin)


###################################################
### code chunk number 9: Vignette.Rnw:220-221 (eval = FALSE)
###################################################
## plot(bin, twoDimensional=TRUE)


###################################################
### code chunk number 10: Vignette.Rnw:225-231
###################################################
pdf("plot_oneD.pdf")
plot(bin)
dev.off()
pdf("plot_twoD.pdf")
plot(bin, twoDimensional=TRUE)
dev.off()


###################################################
### code chunk number 11: Vignette.Rnw:247-249
###################################################
label <- c(rep(0,5), rep(1,5))
bin <- binarize.kMeans(binarizationExample[10,])


###################################################
### code chunk number 12: Vignette.Rnw:254-257 (eval = FALSE)
###################################################
## plot(bin, twoDimensional=TRUE, 
##      col=label+1, pch=label, 
##      showLegend=FALSE)


###################################################
### code chunk number 13: Vignette.Rnw:259-262
###################################################
pdf("plot_bin_with_label.pdf")
plot(bin, twoDimensional=TRUE, col=label+1, pch=label, showLegend=FALSE)
dev.off()


###################################################
### code chunk number 14: Vignette.Rnw:275-277
###################################################
binMatrix <- binarizeMatrix(binarizationExample, 
                            method="kMeans")


###################################################
### code chunk number 15: Vignette.Rnw:281-284
###################################################
binMatrixFDR <- binarizeMatrix(binarizationExample, 
                               method="kMeans",
                               adjustment="fdr")


###################################################
### code chunk number 16: Vignette.Rnw:291-293
###################################################
bin <- binarize.BASC(binarizationExample[1,], method="A")
print(bin)


###################################################
### code chunk number 17: Vignette.Rnw:307-308
###################################################
print(bin@intermediateStrongestSteps)


###################################################
### code chunk number 18: Vignette.Rnw:314-320
###################################################
pdf("stepsA.pdf")
plotStepFunctions(bin, connected=TRUE)
dev.off()
pdf("stepsB.pdf")
plotStepFunctions(binarize.BASC(binarizationExample[1,], method="B"), connected=TRUE)
dev.off()


###################################################
### code chunk number 19: Vignette.Rnw:322-323 (eval = FALSE)
###################################################
## plotStepFunctions(bin)


###################################################
### code chunk number 20: Vignette.Rnw:350-355
###################################################
binMatrix <- binarizeMatrix(binarizationExample, 
                            method="kMeans",
                            adjustment="fdr")
significantRows <- sum(binMatrix[,12] < 0.05)
print(significantRows)


###################################################
### code chunk number 21: Vignette.Rnw:363-369
###################################################
binarizations <- apply(binarizationExample, 1, binarize.BASC, method="A")
pVals <- p.adjust(sapply(binarizations, function(x)
    {
      return(x@p.value)
    }), method="fdr")
  significantRows <- sum(pVals < 0.05)


###################################################
### code chunk number 22: Vignette.Rnw:371-372
###################################################
print(significantRows)


###################################################
### code chunk number 23: Vignette.Rnw:377-383
###################################################
binarizations <- apply(binarizationExample, 1, binarize.BASC, method="B")
pVals <- p.adjust(sapply(binarizations, function(x)
    {
      return(x@p.value)
    }), method="fdr")
significantRows <- sum(pVals < 0.05)


###################################################
### code chunk number 24: Vignette.Rnw:385-386
###################################################
print(significantRows)


###################################################
### code chunk number 25: Vignette.Rnw:405-407
###################################################
tauValues <- seq(0,0.25, 0.05)
print(tauValues)


###################################################
### code chunk number 26: Vignette.Rnw:412-425
###################################################
significantFeatures <- sapply(tauValues, function(tau)
{
  binMatrix <- binarizeMatrix(binarizationExample, 
                              method="BASCB", 
                              adjustment="fdr",
                              tau=tau)
  
  significantRows <- sum(binMatrix[,12] < 0.05)  
  return(significantRows)})

names(significantFeatures) <- tauValues

print(significantFeatures)


