# TODO: Determine how to extend correlation view to multiple assets
# TODO: Create plot distributions function

# This example script compares the numerical and the analytical solution of
# entropy-pooling, see "A. Meucci - Fully Flexible Views: Theory and Practice"
# and associated script S_Main
# Example compares analytical vs. numerical approach to entropy pooling

# Code by A. Meucci, September 2008
# Last version available at www.symmys.com > Teaching > MATLAB

################################################################################
# prior
################################################################################
library(MASS)
# analytical representation
N <- 2 # market dimension (2 assets)
Mu <- zeros(N, 1)
r <- .6
# nxn correlation matrix with correlation 'r' in off-diagonals
Sigma <- (1 - r) * eye(N) + r * ones(N, N)

# numerical representation
J <- 100000 # number of scenarios
p <- ones(J, 1) / J
# distribution centered on (0,0) with variance Sigma
dd <- mvrnorm(J / 2, zeros(N, 1), Sigma)
X <- ones(J, 1) %*% t(Mu) + rbind(dd, -dd) # JxN matrix of scenarios

################################################################################
# views
################################################################################

# location
# long the first and asset and short the second asset produces an expectation
# (of Mu_Q calculated below)
Q <- matrix(c(1, -1), nrow = 1)
Mu_Q <- .5

# scatter
G <- matrix(c(-1, 1), nrow = 1)
Sigma_G <- .5 ^ 2

################################################################################
#  posterior 
################################################################################

# analytical posterior
RevisedMuSigma <- Prior2Posterior(Mu, Q, Mu_Q, Sigma, G, Sigma_G)
Mu_ <- RevisedMuSigma$M_
Sigma_ <- RevisedMuSigma$S_

# numerical posterior
Aeq <- ones(1, J)  # constrain probabilities to sum to one...
beq <- 1

# create views
QX <- X %*% t(Q) # a Jx1 matrix

Aeq <- rbind(Aeq, t(QX))   # ...constrain the first moments...
# QX is a linear combination of vector Q and the scenarios X

beq <- rbind(beq, Mu_Q)

# ...constrain the second moments... 
SecMom <- G %*% Mu_ %*% t(Mu_) %*% t(G) + Sigma_G
# We use Mu_ from analytical result. We do not use Revised Sigma because we are
# testing whether the numerical approach for handling expectations of covariance
# matches the analytical approach.
# TODO: Can we perform this procedure without relying on Mu_ from analytical res
GX <- X %*% t(G)

for (k in 1:nrow(G)) {
  for (l in k:nrow(G)) {
    Aeq <- rbind(Aeq, t(GX[, k] * GX[, l]))
    beq <- rbind(beq, SecMom[k, l])
  }
}

emptyMatrix <- matrix(, nrow = 0, ncol = 0)
# ...compute posterior probabilities
p_ <- EntropyProg(p, emptyMatrix, emptyMatrix, Aeq, beq)$p_

################################################################################
# plots
################################################################################
PlotDistributions(X, p, Mu, Sigma, p_, Mu_, Sigma_)
