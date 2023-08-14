#' @title tdGBH
#' @param p.mat a matrix of p values. Rows: features, columns: outcomes.
#' @param pi0.method a character string of 'lsl','tst' or 'storey'. Default is storey. The name of any method used in estimate.pi0.
#' @param global.pi0.method a character string of 'lsl','tst' or 'storey'. Default is storey. The name of any method used in estimate.pi0.
#' @param shrink a number between 0-1. Default is 0.1. A shrinkage factor to reduce the effects of sampling variation,which is the weighted average of the global and group-specific proportions of null hypotheses.
#' @return The adjusted FDR.
#' @rdname tdGBH
#' @import structSSI
#' @import stats
#' @export tdGBH


tdGBH <- function(p.mat, pi0.method = 'storey',global.pi0.method = 'storey', shrink = 0.1){

  pi0.method <- match.arg(pi0.method)
  global.pi0.method  <- match.arg( global.pi0.method)

  pi0 <- structSSI::estimate.pi0(as.vector(p.mat), method =  global.pi0.method)

  pi0.o <- apply(p.mat, 2, function(x) structSSI::estimate.pi0(x, method = pi0.method))
  pi0.g <- apply(p.mat, 1, function(x) structSSI::estimate.pi0(x, method = pi0.method))

  pi0.o <- (1 - shrink) * pi0.o + shrink * pi0
  pi0.g <- (1 - shrink)  * pi0.g + shrink * pi0

  pi0.o.mat <- t(matrix(pi0.o, nrow = length(pi0.o), ncol = length(pi0.g)))
  pi0.g.mat <- matrix(pi0.g, nrow = length(pi0.g), ncol = length(pi0.o))

  sd.r <- (sd(pi0.o) / sqrt(length(pi0.o)))  / (sd(pi0.o) / sqrt(length(pi0.o)) + sd(pi0.g) / sqrt(length(pi0.g)))
  pi0.og.mat <- sqrt((pi0.o.mat^(2 * sd.r)) * (pi0.g.mat^(2 * (1 - sd.r))))

  pi0 <- mean(pi0.og.mat)
  ws.og.mat <- (1 - pi0.og.mat) / pi0.og.mat
  p.ws.mat <- p.mat / ws.og.mat * (1 - pi0)
  p.adj <- matrix(p.adjust(as.vector(p.ws.mat), 'BH'), length(pi0.g), length(pi0.o))

  return(p.adj)
}

