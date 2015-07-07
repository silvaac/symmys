# This script estimates the prior of a hedge fund return and processes extreme
# views on CVaR according to Entropy Pooling
#
# This script complements the article "Fully Flexible Extreme Views" by
# A. Meucci, D. Ardia, S. Keel available at www.ssrn.com
#	
# The most recent version of this code is available at MATLAB Central -
# File Exchange

# IMPORTANT - This script is about the methodology, not the input data, 
# which has been modified

data(pseudodata)
data(ghqx)

xi <- as.matrix(100 * pseudodata[, 2])
n <- nrow(xi)

# bandwidth
bw <- kernelbw(xi)

# weights
lambda <- log(2) / (n / 2)
wi <- as.matrix(exp(-lambda * (n - seq(from = n, to = 1))))
wi <- as.matrix(apply(wi, 2, rev) / sum(wi))
                
# Prior market model
# kernel density
market <- list()
market$mu  <- mean(xi)
market$pdf <- function (x) kernelpdf(x, xi, bw, wi)
market$cdf <- function (x) kernelcdf(x, xi, bw, wi)
market$inv <- function (x) kernelinv(x, xi, bw, wi)
market$VaR95 <- market$inv(c(0.05))
market$CVaR95 <- integrate(function(x) x * market$pdf(x), -100,
                                                 market$VaR95)$val / 0.05

# numerical (Gauss-Hermite grid) prior

# rescale GH zeros so they belong to [0,1]
tmp <- (ghqx - min(ghqx)) / (max(ghqx) - min(ghqx))
epsilon <- 1e-10
Lower <- market$inv(epsilon)
Upper <- market$inv(1 - epsilon)
X  <- Lower + tmp * (Upper - Lower) # rescale mesh
                
p <- integrateSubIntervals(X, market$cdf)
p <- normalizeProb(p)
J <- nrow(X)
                
# Entropy posterior from extreme view on mean and CVaR
view <- list()
view$mu   <- mean(xi) - 1.0
view$CVaR95 <- market$CVaR95 - 1.0
                
# Netwton Raphson
emptyMatrix <- matrix(,nrow = 0, ncol = 0)
s <- emptyMatrix
idx <- as.matrix(cumsum(p) <= 0.05)
s[1] <- sum(idx)
posteriorEntropy <- EntropyProg(p, t(idx), as.matrix(0.05), rbind(rep(1, J),
                               t(X), t(idx * X)), rbind(1, view$mu,
                               0.05 * view$CVaR95))
KLdiv <- as.matrix(posteriorEntropy$optimizationPerformance$ml)
p_ <- posteriorEntropy$p_

doStop <- 0
i <- 1
while (!doStop) {
  i <- i + 1

  idx <- cbind(matrix(1, 1, s[i - 1]), matrix(0, 1, J - s[i - 1]))
  posteriorEntropy1 <- EntropyProg(p, idx, as.matrix(0.05),
                                  rbind(matrix(1, 1, J), t(X), t(t(idx) * X)),
                                  rbind(1, view$mu, 0.05 * view$CVaR95))
  
  idx <- cbind(matrix(1, 1, s[i - 1] + 1), matrix(0, 1, J - s[i - 1] - 1))
  posteriorEntropy2 <- EntropyProg(p, idx, as.matrix(0.05),
                                  rbind(matrix(1, 1, J), t(X), t(t(idx) * X)),
                                  rbind(1, view$mu, 0.05 * view$CVaR95))
  
  idx <- cbind(matrix(1, 1, s[i - 1] + 2), matrix(0, 1, J - s[i - 1] - 2))
  posteriorEntropy3 <- EntropyProg(p, idx, as.matrix(0.05),
                                  rbind(matrix(1, 1, J), t(X), t(t(idx) * X)),
                                  rbind(1, view$mu, 0.05 * view$CVaR95))
  
  # first difference
  DE  <- posteriorEntropy2$optimizationPerformance$ml -
        posteriorEntropy1$optimizationPerformance$ml

  # second difference
  D2E <- posteriorEntropy3$optimizationPerformance$ml -
        2 * posteriorEntropy2$optimizationPerformance$ml +
        posteriorEntropy1$optimizationPerformance$ml

  # optimal s
  s <- rbind(s, round(s[i - 1] - (DE / D2E)))

  tmp <- emptyMatrix
  idx  <- cbind(matrix(1, 1, s[i]), matrix(0, 1, J - s[i]))
  tempEntropy <- EntropyProg(p, idx, as.matrix(0.05), rbind(matrix(1, 1, J),
                            t(X), t(t(idx) * X)),
                            rbind(1, view$mu, 0.05 * view$CVaR95))
  p_ <- cbind(p_, tempEntropy$p_)
  KLdiv <- rbind(KLdiv, tempEntropy$optimizationPerformance$ml)

  # if change in KLdiv less than one percent, stop
  if(abs((KLdiv[i] - KLdiv[i - 1]) / KLdiv[i - 1])  < 0.01)
    doStop <- 1
}

plot(t(X), p, type = "l", col = "blue", lwd = 1.5, xlab="Returns %",
     ylab = "")
lines(t(X), p_[, ncol(p_)], type = "l", col = "red", lwd = 1.5)
lines(market$mu, 0.0,  type = "p", pch = 16, col = "blue")
lines(market$CVaR95, 0.0,  type = "p", pch = 17, col = "blue")
lines(view$mu, 0.0,  type = "p", pch = 16, col = "red")
lines(view$CVaR95, 0.0,  type = "p", pch = 17, col = "red")
legend("topright", 1.9, c("Prior", "Posterior"),
       col = c("blue", "red"), lty = 1)

