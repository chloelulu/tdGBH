#' @title Two-Dimensional Group Benjamini-Hochberg (tdGBH) Procedure
#'
#' @description Two-Dimensional Group Benjamini-Hochberg (tdGBH) Procedure is a multiple testing procedure that offers enhanced power over traditional FDR control methods when the data presents two-way grouping structure. It reweights the p-values based on the informativeness of the respective grouping directions. It requires a matrix of p-values as input with rows and columns corresponding to the two grouping directions (e.g. genes by cell types).
#'
#' @param p.mat A matrix representing p-values with rows and columns corresponding to the two grouping directions (e.g. genes by cell types).
#'
#' @param pi0.method A character string of 'lsl', 'tst', or 'storey', specifying the method for estimating the group-specific null proportions. The default is 'storey'.
#'
#' @param global.pi0.method A character string of 'lsl', 'tst', or 'storey', specifying the method for estimating the global null proportion. The default is 'storey'.
#'
#' @param shrinkage.method A character string of 'linear' and 'power', indicating whether the shrinkage is on the  linear or power scale. 
#' 
#' @param shrinkage.factor A numeric value between 0 and 1. The factor determines the degree of shrinkage to the global null proportion estimate (0: no shrinkage, 1: maximum shrinkage, i.e. using the global estimate). The default is 0.1.
#'
#' @return Returns the FDR-adjusted p-values.
#'
#' @import stats
#' @author Jun Chen and Lu Yang
#' @references Lu Yang, Pei Wang, Jun Chen. 2dGBH: Two-dimensional Group Benjamini-Hochberg Procedure for False Discovery Rate Control in Two-Way Multiple Testing of Genomic Data
#' @examples
#' data(P)
#' p.adj <- tdGBH(P)
#' @export tdGBH


tdGBH <- function(p.mat, pi0.method = c('storey', 'lsl', 'tst'), global.pi0.method = c('storey', 'lsl', 'tst'), shrinkage.method = c('linear', 'power'), shrinkage.factor = 0.1){

  pi0.method <- match.arg(pi0.method, c( 'storey','lsl','tst'))
  global.pi0.method  <- match.arg(global.pi0.method, c( 'storey','lsl','tst'))
  shrinkage.method <- match.arg(shrinkage.method, c('linear', 'power'))
  
  pi0 <- estimate.pi0(as.vector(p.mat), method =  global.pi0.method)

  pi0.o <- apply(p.mat, 2, function(x) estimate.pi0(x, method = pi0.method))
  pi0.g <- apply(p.mat, 1, function(x) estimate.pi0(x, method = pi0.method))
  
  if(shrinkage.method == 'linear'){
    pi0.o <- (1 - shrinkage.factor) * pi0.o + shrinkage.factor * pi0
    pi0.g <- (1 - shrinkage.factor)  * pi0.g + shrinkage.factor * pi0
  }
  
  if(shrinkage.method == 'power'){
    pi0.o <-  (pi0.o^(1 - shrinkage.factor)) * (pi0^shrinkage.factor)
    pi0.g <- (pi0.g^(1 - shrinkage.factor)) * (pi0^shrinkage.factor)
  }

  pi0.o.mat <- t(matrix(pi0.o, nrow = length(pi0.o), ncol = length(pi0.g)))
  pi0.g.mat <- matrix(pi0.g, nrow = length(pi0.g), ncol = length(pi0.o))

  sd.r <- (sd(pi0.o) / sqrt(length(pi0.o)))  / (sd(pi0.o) / sqrt(length(pi0.o)) + sd(pi0.g) / sqrt(length(pi0.g)))
  pi0.og.mat <- sqrt((pi0.o.mat^(2 * sd.r)) * (pi0.g.mat^(2 * (1 - sd.r))))

  pi0 <- mean(pi0.og.mat)
  ws.og.mat <- (1 - pi0.og.mat) / pi0.og.mat
  p.ws.mat <- p.mat / ws.og.mat * (1 - pi0)

  p.adj <- matrix(p.adjust(as.vector(p.ws.mat), 'BH'), length(pi0.g), length(pi0.o), dimnames = list(rownames(p.mat),colnames(p.mat)))
  return(p.adj)
}


pi0.tail.p <- function(lambda, p.values){
  num <- length(which(p.values >= lambda))
  denom <- length(p.values)*(1 - lambda)
  min(num/denom, 1)
}


pi0.tst <- function(p.val, alpha = 0.05){
  alpha.prime <- alpha/(1 + alpha)
  n_g <- length(p.val)
  adjustment <- mt.rawp2adjp(p.val, proc = "BH")
  rejected <- mt.reject(adjustment$adjp, alpha.prime)
  n.rejected <- rejected$r[,2]
  (n_g - n.rejected) / n_g
}

pi0.lsl <- function(p.val){
  p.val <- sort(p.val)
  n_g <- length(p.val)
  
  i <- 1
  while(TRUE){
    if(i >= 2){
      l_g.i.prev <- l_g.i
    } else {
      l_g.i.prev <- Inf
    }
    if(p.val[i] < 1){
      l_g.i <- (n_g + 1 - i)/(1 - p.val[i])
    } else {
      return(l_g.i.prev) 
    }
    if(l_g.i > l_g.i.prev || i == length(p.val)){
      pi0 <- (floor(l_g.i) + 1)/n_g
      pi0 <- min(pi0, 1)
      return(pi0)
    }
    i <- i + 1
  }
}

estimate.pi0 <- function (pvalues, method, alpha = 0.05, lambda = 0.5) {
  method <- tolower(method)
  matched.method <- match.arg(method, c("tst", "lsl", "storey"))
  if (matched.method == "tst") {
    return(pi0.tst(pvalues, alpha))
  }
  else if (matched.method == "lsl") {
    return(pi0.lsl(pvalues))
  }
  else if (matched.method == "storey") {
    return(pi0.tail.p(lambda, pvalues))
  }
}



mt.rawp2adjp <- function (rawp, proc = c("Bonferroni", "Holm", "Hochberg", "SidakSS", 
                         "SidakSD", "BH", "BY", "ABH", "TSBH"), alpha = 0.05, na.rm = FALSE) 
{
  m <- length(rawp)
  if (na.rm) {
    mgood <- sum(!is.na(rawp))
  }
  else {
    mgood <- m
  }
  n <- length(proc)
  a <- length(alpha)
  index <- order(rawp)
  h0.ABH <- NULL
  h0.TSBH <- NULL
  spval <- rawp[index]
  adjp <- matrix(0, m, n + 1)
  dimnames(adjp) <- list(NULL, c("rawp", proc))
  adjp[, 1] <- spval
  if (is.element("TSBH", proc)) {
    TS.spot <- which(proc == "TSBH")
    TSBHs <- paste("TSBH", alpha, sep = "_")
    newprocs <- append(proc, TSBHs, after = TS.spot)
    newprocs <- newprocs[newprocs != "TSBH"]
    adjp <- matrix(0, m, n + a)
    dimnames(adjp) <- list(NULL, c("rawp", newprocs))
    adjp[, 1] <- spval
    tmp <- spval
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood/i) * spval[i], 
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i])) 
        tmp[i] <- NA
    }
    h0.TSBH <- rep(0, length(alpha))
    names(h0.TSBH) <- paste("h0.TSBH", alpha, sep = "_")
    for (i in 1:length(alpha)) {
      h0.TSBH[i] <- mgood - sum(tmp < alpha[i]/(1 + alpha[i]), 
                                na.rm = TRUE)
      adjp[, TS.spot + i] <- tmp * h0.TSBH[i]/mgood
    }
  }
  if (is.element("Bonferroni", proc)) {
    tmp <- mgood * spval
    tmp[tmp > 1] <- 1
    adjp[, "Bonferroni"] <- tmp
  }
  if (is.element("Holm", proc)) {
    tmp <- spval
    tmp[1] <- min(mgood * spval[1], 1)
    for (i in 2:m) tmp[i] <- max(tmp[i - 1], min((mgood - 
                                                    i + 1) * spval[i], 1))
    adjp[, "Holm"] <- tmp
  }
  if (is.element("Hochberg", proc)) {
    tmp <- spval
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood - i + 1) * spval[i], 
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i])) 
        tmp[i] <- NA
    }
    adjp[, "Hochberg"] <- tmp
  }
  if (is.element("SidakSS", proc)) 
    adjp[, "SidakSS"] <- 1 - (1 - spval)^mgood
  if (is.element("SidakSD", proc)) {
    tmp <- spval
    tmp[1] <- 1 - (1 - spval[1])^mgood
    for (i in 2:m) tmp[i] <- max(tmp[i - 1], 1 - (1 - spval[i])^(mgood - 
                                                                   i + 1))
    adjp[, "SidakSD"] <- tmp
  }
  if (is.element("BH", proc)) {
    tmp <- spval
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood/i) * spval[i], 
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i])) 
        tmp[i] <- NA
    }
    adjp[, "BH"] <- tmp
  }
  if (is.element("BY", proc)) {
    tmp <- spval
    a <- sum(1/(1:mgood))
    tmp[m] <- min(a * spval[m], 1)
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood * a/i) * spval[i], 
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i])) 
        tmp[i] <- NA
    }
    adjp[, "BY"] <- tmp
  }
  if (is.element("ABH", proc)) {
    tmp <- spval
    h0.m <- rep(0, mgood)
    for (k in 1:mgood) {
      h0.m[k] <- (mgood + 1 - k)/(1 - spval[k])
    }
    grab <- min(which(diff(h0.m, na.rm = TRUE) > 0), na.rm = TRUE)
    h0.ABH <- ceiling(min(h0.m[grab], mgood))
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood/i) * spval[i], 
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i])) 
        tmp[i] <- NA
    }
    adjp[, "ABH"] <- tmp * h0.ABH/mgood
  }
  list(adjp = adjp, index = index, h0.ABH = h0.ABH[1], h0.TSBH = h0.TSBH[1:length(alpha)])
}



mt.reject <- function (adjp, alpha) 
{
  which <- adjp <= alpha[1]
  dimnames(which) <- dimnames(adjp)
  if (is.matrix(adjp)) {
    r <- matrix(0, length(alpha), ncol(adjp))
    for (i in 1:length(alpha)) r[i, ] <- colSums(adjp <= 
                                                   alpha[i])
    dimnames(r) <- list(alpha, dimnames(adjp)[[2]])
  }
  if (!is.matrix(adjp)) {
    r <- rep(0, length(alpha))
    for (i in 1:length(alpha)) r[i] <- sum(adjp <= alpha[i])
  }
  list(r = r, which = which)
}

