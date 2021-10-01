#!/usr/bin/env Rscript

## Modified version of perm.t.test from GmAMisc
## Original version: https://CRAN.R-project.org/package=GmAMisc
## Modified by: Jana Biermann, PhD

perm.t.test.mod <- function(data, format, sample1.lab = NULL, sample2.lab = NULL, B = 999, pathway,plot) {
  #options(scipen = 999)
  if (format == "long") {
    unstacked.data <- utils::unstack(data)
    sample1 <- unstacked.data[[1]]
    sample2 <- unstacked.data[[2]]
  }
  else {
    sample1 <- data[, 1]
    sample2 <- data[, 2]
  }
  if (is.null(sample1.lab) == TRUE) {
    sample1.lab <- "smpl 1"
    sample2.lab <- "smpl 2"
  }
  else {
  }
  n1 <- length(sample1)
  n2 <- length(sample2)
  mean1 <- mean(sample1)
  mean2 <- mean(sample2)
  error1 <- qnorm(0.975) * sd(sample1)/sqrt(n1)
  error2 <- qnorm(0.975) * sd(sample2)/sqrt(n2)
  sample1_lci <- mean1 - error1
  sample1_uci <- mean1 + error1
  sample2_lci <- mean2 - error2
  sample2_uci <- mean2 + error2
  p.equal.var <- t.test(sample1, sample2, var.equal = TRUE)$p.value
  p.unequal.var <- t.test(sample1, sample2, var.equal = FALSE)$p.value
  pooledData <- c(sample1, sample2)
  size.sample1 <- length(sample1)
  size.sample2 <- length(sample2)
  size.pooled <- size.sample1 + size.sample2
  nIter <- B
  meanDiff <- numeric(nIter + 1)
  meanDiff[1] <- mean1 - mean2
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (i in 2:B) {
    index <- sample(1:size.pooled, size = size.sample1, replace = F)
    sample1.perm <- pooledData[index]
    sample2.perm <- pooledData[-index]
    meanDiff[i] <- mean(sample1.perm) - mean(sample2.perm)
    setTxtProgressBar(pb, i)
  }
  p.lowertail <- (1 + sum(meanDiff[-1] < meanDiff[1]))/(1 + B)
  p.uppertail <- (1 + sum(meanDiff[-1] > meanDiff[1]))/(1 + B)
  two.sided.p <- 2 * min(p.lowertail, p.uppertail)
  p.wilcox <- wilcox.test(sample1, sample2)$p.value
  
  if(plot==T){
    graphics::hist(meanDiff, main = paste0(pathway, "\nDistribution of permuted mean differences\n(number of permutations: ", B, "; q10 and q90 excluded)"), 
                   xlab = "", sub = paste0(sample1.lab, " (n: ", n1, ") (95% CI lower bound., mean, 95% CI upper bound.): ", round(sample1_lci,2), ", ", round(mean1, 2), ", ", round(sample1_uci, 2), "\n", 
                                           sample2.lab, " (n: ", n2, ") (95% CI lower bound., mean, 95% CI upper bound.): ", round(sample2_lci,2), ", ", round(mean2, 2), ", ", round(sample2_uci, 2), "\nobs. mean diff. (solid dot): ",
                                           round(meanDiff[1], 2), "; perm. P mean ", sample1.lab, " < ", sample2.lab, ": ", round(p.lowertail, 4), 
                                           "; perm. P mean ", sample1.lab, " > ", sample2.lab, ": ", round(p.uppertail, 4), 
                                           "\nperm. P (2-sided): ", round(two.sided.p, 4), 
                                           "; regular t-test P (2-sided): ", round(p.equal.var, 4), " (equal var); ", round(p.unequal.var, 4), " (unequal var)"), cex.main = 0.85, cex.sub = 0.7)
    rug(meanDiff, col = "#0000FF")
    abline(v = stats::quantile(meanDiff, 0.025), lty = 2, col = "blue")
    abline(v = stats::quantile(meanDiff, 0.975), lty = 2, col = "blue")
    points(x = meanDiff[1], y = 0, pch = 20, col = "black")
    points(x = mean(meanDiff[-1]), y = 0, pch = 1, col = "black")
  }
  
  res_perm<<-data.frame(pathway, sample1.lab, n1, sample1_lci, mean1, sample1_uci, sample2.lab, n2, 
                        sample2_lci, mean2, sample2_uci, meanDiff[1], p.lowertail, p.uppertail,
                        two.sided.p, p.equal.var, p.unequal.var, p.wilcox)
  return(res_perm)
}