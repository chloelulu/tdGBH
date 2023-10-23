#' @title Two-Dimensional Group Benjamini-Hochberg (tdGBH) Procedure
#'
#' @description This function applies the two-dimensional group Benjamini-Hochberg (tdGBH) procedure to control the false discovery rate.
#'
#' @param p.mat A matrix representing p-values where rows correspond to features and columns to outcomes.
#'
#' @param pi0.method A character string specifying the method for \link[structSSI]{estimate.pi0}, chosen from 'lsl', 'tst', or 'storey'. The default is 'storey'.
#'
#' @param global.pi0.method A character string that indicates the method used for the global \link[structSSI]{estimate.pi0}. Possible values are 'lsl', 'tst', or 'storey'. The default is 'storey'.
#'
#' @param weight.method A character string specifying the method for add the weight, chosen from 'geo', 'ari', or 'new', represents geometric mean,  arithmetic mean, and tdgbh method. The default is 'new'.
#' 
#' @param shrink A numeric value between 0 and 1, serving as a shrinkage factor. It is employed to mitigate the impact of sampling variability. The factor determines the weighted average between global and group-specific proportions of null hypotheses. The default is 0.1.
#'
#' @param renorm A logical value indicating whether renormalization should be performed on the weights, which ensures FDR control for non data-driven weights. The default is FALSE.
#' 
#' @return Returns the adjusted p-values.
#'
#' @import stats
#' @author Lu Yang, Pei Wang and Jun Chen
#' @references Lu Yang, Pei Wang, Jun Chen. 2dGBH: Two-dimensional Group Benjamini-Hochberg Procedure for False Discovery Rate Control in Two-Way Multiple Testing.
#' @examples
#' data(P)
#' p.adj <- tdGBH(P)
#' @export tdGBH


tdGBH <- function(p.mat, pi0.method = 'storey', global.pi0.method = 'storey', weight.method = 'new', shrink = 0.1, renorm=F){

  pi0.method <- match.arg(pi0.method)
  global.pi0.method  <- match.arg( global.pi0.method)

  pi0 <- estimate.pi0(as.vector(p.mat), method =  global.pi0.method)

  pi0.o <- apply(p.mat, 2, function(x) estimate.pi0(x, method = pi0.method))
  pi0.g <- apply(p.mat, 1, function(x) estimate.pi0(x, method = pi0.method))

  pi0.o <- (1 - shrink) * pi0.o + shrink * pi0
  pi0.g <- (1 - shrink)  * pi0.g + shrink * pi0

  pi0.o.mat <- t(matrix(pi0.o, nrow = length(pi0.o), ncol = length(pi0.g)))
  pi0.g.mat <- matrix(pi0.g, nrow = length(pi0.g), ncol = length(pi0.o))

  if(weight.method == 'geo') pi0.og.mat <- sqrt(pi0.o.mat * pi0.g.mat)
  if(weight.method == 'ari') pi0.og.mat <- (pi0.o.mat + pi0.g.mat) / 2
  if(weight.method == 'new') {
    sd.r <- (sd(pi0.o) / sqrt(length(pi0.o)))  / (sd(pi0.o) / sqrt(length(pi0.o)) + sd(pi0.g) / sqrt(length(pi0.g)))
    pi0.og.mat <- sqrt((pi0.o.mat^(2 * sd.r)) * (pi0.g.mat^(2 * (1 - sd.r))))
  }
  
  pi0 <- mean(pi0.og.mat)
  ws.og.mat <- (1 - pi0.og.mat) / pi0.og.mat
  if (renorm == TRUE) {
    ws <- ws.og.mat / (1 - pi0)
    ws <- length(p.mat) / sum(ws) * ws
    p.ws.mat <- p.mat / ws
  } else {
    p.ws.mat <- p.mat / ws.og.mat * (1 - pi0)
  }
  
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

