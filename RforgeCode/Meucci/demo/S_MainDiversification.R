#' @title computes Diversification Distribution and the Effective Number of Bets
#'
#' @description  computes Diversification Distribution and the Effective Number 
#' of Bets with 1) Principal Components Bets, 2) Minimum Torsion Bets,and 
#' computes the Marginal Contribution to Risk for an equal-weight portfolio of 
#' stocks in the S&P500, as described in A. Meucci, A. Santangelo, R. Deguest - 
#' "Measuring Portfolio Diversification Based on Optimized Uncorrelated Factors"
#'
#' @references
#' A. Meucci, A. Santangelo, R. Deguest - "Measuring Portfolio Diversification 
#' Based on Optimized Uncorrelated Factors" \url{http://symmys.com/node/599}
#'
#' See Meucci's script "S_MainDiversification.m"
#
#' @author Xavier Valls \email{xaviervallspla@@gmail.com}

data(linRet)

# equally weighted exposures (weights) to factors (returns)
n_ <- nrow(linRet);
b <- ones( n_, 1 ) / n_;

#Sample covariance matrix
Sigma <- cov(t(linRet)) ;

# PCA decomposition
e <- eigen(Sigma);

# PCA torsion matrix and exposures for ew portfolio
t_PC <- Torsion(Sigma, 'pca');

# Minimum-Torsion matrix and exposures for ew portfolio
t_MT <- Torsion(Sigma, 'minimum-torsion', 'exact');

# Diversification Distribition and NEB using PCA torsion matrix
ENB <- EffectiveBets(b, Sigma, t_PC);
ENB_PC <- ENB$enb
DiverDistr_PC <- ENB$p
# Diversification Distribition and NEB using Minimum-Torsion matrix
ENB <- EffectiveBets(b, Sigma, t_MT);
ENB_MT <- ENB$enb
DiverDistr_MT <- ENB$p
# Marginal Risk Contribution
MargContrs <- b * (Sigma %*% b)/(t(b) %*% Sigma %*% b)[1];
