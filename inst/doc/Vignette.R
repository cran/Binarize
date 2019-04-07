### R code from vignette source 'Vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Vignette.Rnw:36-37
###################################################
set.seed(13579)


###################################################
### code chunk number 2: Vignette.Rnw:195-196 (eval = FALSE)
###################################################
## install.packages("Binarize")


###################################################
### code chunk number 3: Vignette.Rnw:199-200
###################################################
library("Binarize")


###################################################
### code chunk number 4: Vignette.Rnw:213-214
###################################################
data(binarizationExample)


###################################################
### code chunk number 5: Vignette.Rnw:220-227
###################################################
pdf("density.pdf")
par(mar=c(2,2,1,1))
#plot(density(binarizationExample[1,]),main="")
#abline(v=mean(binarizationExample[1,]), lty="dashed")
plot(function(x)dnorm(x,mean=0,sd=1)+dnorm(x,mean=10,sd=1),xlim=c(-5,15),main="")
abline(v=5, lty="dashed")
dev.off()


###################################################
### code chunk number 6: Vignette.Rnw:241-243
###################################################
bin <- binarize.kMeans(binarizationExample[1,])
print(bin)


###################################################
### code chunk number 7: Vignette.Rnw:250-251
###################################################
print(bin@binarizedMeasurements)


###################################################
### code chunk number 8: Vignette.Rnw:255-256 (eval = FALSE)
###################################################
## plot(bin)


###################################################
### code chunk number 9: Vignette.Rnw:259-260 (eval = FALSE)
###################################################
## plot(bin, twoDimensional=TRUE)


###################################################
### code chunk number 10: Vignette.Rnw:264-270
###################################################
pdf("plot_oneD.pdf")
plot(bin)
dev.off()
pdf("plot_twoD.pdf")
plot(bin, twoDimensional=TRUE)
dev.off()


###################################################
### code chunk number 11: Vignette.Rnw:286-288
###################################################
label <- c(rep(0,5), rep(1,5))
bin <- binarize.kMeans(binarizationExample[10,])


###################################################
### code chunk number 12: Vignette.Rnw:293-296 (eval = FALSE)
###################################################
## plot(bin, twoDimensional=TRUE, 
##      col=label+1, pch=label, 
##      showLegend=FALSE)


###################################################
### code chunk number 13: Vignette.Rnw:298-301
###################################################
pdf("plot_bin_with_label.pdf")
plot(bin, twoDimensional=TRUE, col=label+1, pch=label, showLegend=FALSE)
dev.off()


###################################################
### code chunk number 14: Vignette.Rnw:314-316
###################################################
binMatrix <- binarizeMatrix(binarizationExample, 
                            method="kMeans")


###################################################
### code chunk number 15: Vignette.Rnw:320-323
###################################################
binMatrixFDR <- binarizeMatrix(binarizationExample, 
                               method="kMeans",
                               adjustment="fdr")


###################################################
### code chunk number 16: Vignette.Rnw:330-332
###################################################
bin <- binarize.BASC(binarizationExample[1,], method="A")
print(bin)


###################################################
### code chunk number 17: Vignette.Rnw:346-347
###################################################
print(bin@intermediateStrongestSteps)


###################################################
### code chunk number 18: Vignette.Rnw:353-359
###################################################
pdf("stepsA.pdf")
plotStepFunctions(bin, connected=TRUE)
dev.off()
pdf("stepsB.pdf")
plotStepFunctions(binarize.BASC(binarizationExample[1,], method="B"), connected=TRUE)
dev.off()


###################################################
### code chunk number 19: Vignette.Rnw:361-362 (eval = FALSE)
###################################################
## plotStepFunctions(bin)


###################################################
### code chunk number 20: Vignette.Rnw:386-387
###################################################
data(trinarizationExample)


###################################################
### code chunk number 21: Vignette.Rnw:393-395
###################################################
tri <- TASC(trinarizationExample[1,], method="A")
print(tri)


###################################################
### code chunk number 22: Vignette.Rnw:423-424
###################################################
print(tri@intermediateStrongestSteps)


###################################################
### code chunk number 23: Vignette.Rnw:431-437
###################################################
pdf("triA.pdf")
par(mfrow = c(1,2), mar = c(2,2,1,1))
plotStepFunctions(tri, connected=TRUE)
par(mar = c(2,2,1,1))
plot(tri, twoDimensional = TRUE)
dev.off()


###################################################
### code chunk number 24: Vignette.Rnw:439-441 (eval = FALSE)
###################################################
## plotStepFunctions(tri)
## plot(tri, twoDimensional = TRUE)


###################################################
### code chunk number 25: Vignette.Rnw:469-474
###################################################
binMatrix <- binarizeMatrix(binarizationExample, 
                            method="kMeans",
                            adjustment="fdr")
significantRows <- sum(binMatrix[,12] < 0.05)
print(significantRows)


###################################################
### code chunk number 26: Vignette.Rnw:482-488
###################################################
binarizations <- apply(binarizationExample, 1, binarize.BASC, method="A")
pVals <- p.adjust(sapply(binarizations, function(x)
    {
      return(x@p.value)
    }), method="fdr")
  significantRows <- sum(pVals < 0.05)


###################################################
### code chunk number 27: Vignette.Rnw:490-491
###################################################
print(significantRows)


###################################################
### code chunk number 28: Vignette.Rnw:496-502
###################################################
binarizations <- apply(binarizationExample, 1, binarize.BASC, method="B")
pVals <- p.adjust(sapply(binarizations, function(x)
    {
      return(x@p.value)
    }), method="fdr")
significantRows <- sum(pVals < 0.05)


###################################################
### code chunk number 29: Vignette.Rnw:504-505
###################################################
print(significantRows)


###################################################
### code chunk number 30: Vignette.Rnw:524-526
###################################################
tauValues <- seq(0,0.25, 0.05)
print(tauValues)


###################################################
### code chunk number 31: Vignette.Rnw:531-544
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


