# This script compares the performance of plain Monte Carlo versus grid in 
# applying Entropy Pooling to process extreme views and complements the article
# "Fully Flexible Extreme Views" A. Meucci, D. Ardia, S. Keel available at 
# www.ssrn.com.
#
# The most recent version of this code is available at MATLAB Central - File 
# Exchange

library(matlab)
data("ghqx")

################################################################################
# Prior market model
################################################################################
# analytical (normal) prior
emptyMatrix <- matrix(nrow = 0, ncol = 0)
market <- list()
market$mu   <- 0.0
market$sig2 <- 1.0
market$pdf <- function(x) dnorm(x, mean = market$mu, sd = sqrt(market$sig2))
market$cdf <- function(x) pnorm(x, mean = market$mu, sd = sqrt(market$sig2))
market$rnd <- function(x) rnorm(x, mean = market$mu, sd = sqrt(market$sig2))
market$inv <- function(x) qnorm(x, mean = market$mu, sd = sqrt(market$sig2))

# numerical (Monte Carlo) prior
monteCarlo <- list()
monteCarlo$J <- 100000
monteCarlo$X <- market$rnd(monteCarlo$J)
monteCarlo$p <- normalizeProb(1 / monteCarlo$J * ones(monteCarlo$J, 1))

# numerical (Gauss-Hermite grid) prior
ghqMesh <- list()

# rescale GH zeros so they belong to [0,1]
tmp <- (ghqx - min(ghqx)) / (max(ghqx) - min(ghqx))
epsilon <- 1e-10
Lower <- market$inv(epsilon)
Upper <- market$inv(1 - epsilon)
ghqMesh$X  <- Lower + tmp * (Upper - Lower) # rescale mesh

p <- integrateSubIntervals(ghqMesh$X, market$cdf)
ghqMesh$p <- normalizeProb(p)
ghqMesh$J <- nrow(ghqMesh$X)

################################################################################
# Entropy posterior from extreme view on expectation
################################################################################
# view of the analyst
view <- list()
view$mu <- -3.0

# analytical (known since normal model has analytical solution)
truePosterior <- list()
truePosterior <- Prior2Posterior(market$mu, 1, view$mu, market$sig2, 0)
truePosterior$pdf <- function(x) dnorm(x, truePosterior$M_,
       								   sqrt(truePosterior$S_))

# numerical (Monte Carlo)
Aeq <- rbind(ones(1, monteCarlo$J), t(monteCarlo$X))
beq <- rbind(1, view$mu)
monteCarloOptimResult <- EntropyProg(monteCarlo$p, emptyMatrix, emptyMatrix,
 									 Aeq, beq)

monteCarlo$p_ <- monteCarloOptimResult$p_
monteCarlo$KLdiv <- monteCarloOptimResult$optimizationPerformance$ml

# numerical (Gaussian-Hermite grid)
Aeq <- rbind(ones(1, ghqMesh$J), t(ghqMesh$X))
beq <- rbind(1, view$mu)
ghqMeshOptimResult <- EntropyProg(ghqMesh$p, emptyMatrix, emptyMatrix, Aeq, beq)

ghqMesh$p_ <- ghqMeshOptimResult$p_
ghqMesh$KLdiv <- ghqMeshOptimResult$optimizationPerformance$ml

################################################################################
# Plots
################################################################################
xmin <- min(ghqMesh$X)
xmax <- max(ghqMesh$X)
ymax <- 1.0
xmesh <- t(linspace(xmin, xmax, ghqMesh$J))

# Monte Carlo
dev.new()
plotDataMC <- PHist(monteCarlo$X, monteCarlo$p_, 50, main = "Monte Carlo",
    				xlim = c(xmin, xmax), ylim = c(0, ymax))
lines(xmesh, market$pdf(xmesh), type = "l", col = "blue")
lines(xmesh, truePosterior$pdf(xmesh),  type = "l", col = "red")
lines(0.0, 0.0,  type = "p", pch = 17, col = "blue")
lines(view$mu, 0.0,  type = "p", pch = 17, col = "red")

# Gauss Hermite Grid
dev.new()
plotDataGHQ <- PHist(data.matrix(ghqMesh$X), ghqMesh$p_, 50,
                     main = "Gauss-Hermite grid",
    				 xlim = c(xmin, xmax), ylim = c(0, ymax))
lines(xmesh, market$pdf(xmesh), type = "l", col = "blue")
lines(xmesh, truePosterior$pdf(xmesh),  type = "l", col = "red")
lines(0.0, 0.0,  type = "p", pch = 17, col = "blue")
lines(view$mu, 0.0,  type = "p", pch = 17, col = "red")
