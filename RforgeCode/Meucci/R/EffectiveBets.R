#' @title computes the Effective Number of Bets and the Diversification 
#' distribution
#'
#' @description  computes the Effective Number of Bets and the Diversification 
#' distribution, as described in A. Meucci, A. Santangelo, R. Deguest - 
#' "Measuring Portfolio Diversification Based on Optimized Uncorrelated Factors"
#'
#' @param  b         [vector] (n_ x 1) exposures
#' @param  Sigma     [matrix] (n_ x n_) covariance matrix
#' @param  t         [matrix] (n_ x n_) torsion matrix
#'  
#' @return enb       [scalar] Effetive Number of Bets
#' @return p         [vector] (n_ x 1) diversification distribution
#'
#' @references
#' A. Meucci, A. Santangelo, R. Deguest - "Measuring Portfolio Diversification 
#' Based on Optimized Uncorrelated Factors" \url{http://symmys.com/node/599}
#'
#' See Meucci's script "EffectiveBets.m"
#
#' @author Xavier Valls \email{xaviervallspla@@gmail.com}
#' @export

EffectiveBets <- function(b, Sigma, t) {
  p <- mldivide(t(t), b) * (t %*% Sigma %*% b) / (t(b) %*% Sigma %*% b)[1];
  enb <- exp(-sum(p * log(1 + (p - 1) * (p > 1e-5))));

  return(list(p = p, enb = enb));
}
