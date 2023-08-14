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
#' @param shrink A numeric value between 0 and 1, serving as a shrinkage factor. It is employed to mitigate the impact of sampling variability. The factor determines the weighted average between global and group-specific proportions of null hypotheses. The default is 0.1.
#'
#' @return Returns the adjusted p-values.
#'
#' @import structSSI
#' @import stats
#' @author Lu Yang and Jun Chen
#' @references Lu Yang, Jun Chen. 2dGBH: Two-dimensional Group Benjamini-Hochberg Procedure for False Discovery Rate Control in Two-Way Multiple Testing.
#' @examples
#' data(P)
#' p.adj <- tdGBH(P)
#' @export 



tdGBH <- function(p.mat, pi0.method = 'storey', global.pi0.method = 'storey', shrink = 0.1){

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
  p.adj <- matrix(p.adjust(as.vector(p.ws.mat), 'BH'), length(pi0.g), length(pi0.o), dimnames = list(rownames(p.mat),colnames(p.mat)))
  return(p.adj)
}

