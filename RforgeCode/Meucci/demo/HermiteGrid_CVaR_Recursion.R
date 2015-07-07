# This script illustrates the discrete Newton recursion  to process views on
# CVaR according to Entropy Pooling and complements the article "Fully Flexible
#	Extreme Views" by A. Meucci, D. Ardia, S. Keel available at www.ssrn.com.
# The most recent version of this code is available at
# MATLAB Central - File Exchange

data(gqhx)

# Prior market model (normal) on grid
emptyMatrix <- matrix(nrow = 0, ncol = 0)
market.mu   <- 0.0
market.sig2 <- 1.0
market.pdf <- function(x) dnorm(x, mean = market.mu, sd = sqrt(market.sig2))
market.cdf <- function(x) pnorm(x, mean = market.mu, sd = sqrt(market.sig2))
market.rnd <- function(x) rnorm(x, mean = market.mu, sd = sqrt(market.sig2))
market.inv <- function(x) qnorm(x, mean = market.mu, sd = sqrt(market.sig2))
market.VaR95 <- market.inv(0.05)
market.CVaR95 <- integrate(function(x) (x * market.pdf(x)), -100,
                                       market.VaR95)$val / 0.05
# rescale GH zeros so they belong to [0,1]
tmp <- (ghqx - min(ghqx)) / (max(ghqx) - min(ghqx))
epsilon <- 1e-10
Lower <- market.inv(epsilon)
Upper <- market.inv(1 - epsilon)
X <- Lower + tmp * (Upper - Lower) # rescale mesh

p <- integrateSubIntervals(X, market.cdf)
p <- normalizeProb(p)
J <- nrow(X)

# Entropy posterior from extreme view on CVaR: brute-force approach

# view of the analyst
view.CVaR95 <- -3.0

# Iterate over different VaR95 levels
nVaR95 <- 100
VaR95  <- seq(view.CVaR95, market.VaR95, len=nVaR95)
p_     <- matrix(NaN, nrow = J, ncol = nVaR95)
s_     <- matrix(NaN, nrow = nVaR95, ncol = 1)
KLdiv  <- matrix(NaN, nrow = nVaR95, ncol = 1)

for (i in 1:nVaR95) {
  idx <- as.matrix(X <= VaR95[i])
  s_[i] <- sum(idx)
  posteriorEntropy <- EntropyProg(p, t(idx),  as.matrix(0.05), rbind(rep(1, J),
                                 t(idx * X)), rbind(1, 0.05 * view.CVaR95))
  p_[,i] <- posteriorEntropy$p_
  # when maximizing, objective turns negative. Bringing it back to normal
  KLdiv[i] <- -posteriorEntropy$optimizationPerformance$ml
}

# Display results
dev.new()
par(mfrow = c(2,1))
plot(s_, KLdiv, "l", col = "blue", xlab = "s", ylab = "KL divergence")
dummy <- min(KLdiv)
idxMin <- which.min(KLdiv)
lines(s_[idxMin], KLdiv[idxMin], type = "b")

tmp <- p_[, idxMin]
tmp <- tmp / sum(tmp)
plot(X, tmp, type = "l", col = "blue", lwd = 1.5, xlab="",
     ylab = "", ylim = c(0, 0.006))
x <- seq(min(X), max(X), len = J)
tmp <- market.pdf(x)
tmp <- tmp / sum(tmp)
lines(x, tmp, type = "l", col = "red", lwd = 1.5)
lines(market.CVaR95, 0.0,  type = "p", pch = 17, col = "blue")
lines(view.CVaR95, 0.0,  type = "p", pch = 17, col = "red")
legend("topright", 1.9, c("Prior", "Posterior"),
       col = c("blue", "red"), lty = 1)

# Entropy posterior from extreme view on CVaR: Newton Raphson approach

s <- emptyMatrix

# initial value
idx <- as.matrix(cumsum(p) <= 0.05)
s[1] <- sum(idx)
posteriorEntropy <- EntropyProg(p, t(idx), as.matrix(0.05), rbind(rep(1, J),
                               t(idx * X)), rbind(1, 0.05 * view.CVaR95))
KLdiv <- as.matrix(-posteriorEntropy$optimizationPerformance$ml)
p_ <- posteriorEntropy$p_

# iterate
doStop <- 0
i <- 1
while (!doStop) {
  i <- i + 1

  idx <- cbind(matrix(1, 1, s[i - 1]), matrix(0, 1, J - s[i - 1]))
  posteriorEntropy1 <- EntropyProg(p, idx, as.matrix(0.05),
                                  rbind(matrix(1, 1, J), t(t(idx) * X)),
                                  rbind(1, 0.05 * view.CVaR95))
  # [dummy, KLdiv_s]  = optimizeEntropy(p, [idx'; (idx .* X)'],
  #  [0.05; 0.05 * view.CVaR95], [ones(1, J); X'], [1; view.mu]);

  idx <- cbind(matrix(1, 1, s[i - 1] + 1), matrix(0, 1, J - s[i - 1] - 1))
  posteriorEntropy2 <- EntropyProg(p, idx, as.matrix(0.05),
                                  rbind(matrix(1, 1, J), t(t(idx) * X)),
                                  rbind(1, 0.05 * view.CVaR95))
  # [dummy, KLdiv_s1] = optimizeEntropy(p, [idx'; (idx .* X)'],
  #  [0.05; 0.05 * view.CVaR95], [ones(1, J); X'], [1; view.mu]);

  idx <- cbind(matrix(1, 1, s[i - 1] + 2), matrix(0, 1, J - s[i - 1] - 2))
  posteriorEntropy3 <- EntropyProg(p, idx, as.matrix(0.05),
                                  rbind(matrix(1, 1, J), t(t(idx) * X)),
                                  rbind(1, 0.05 * view.CVaR95))
  # [dummy, KLdiv_s2] = optimizeEntropy(p, [idx'; (idx .* X)'],
  #  0.05; 0.05 * view.CVaR95], [ones(1, J); X'], [1; view.mu]);

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
                            t(t(idx) * X)), rbind(1, 0.05 * view.CVaR95))
  # [tmp.p_, tmp.KLdiv] = optimizeEntropy(p, [idx'; (idx .* X)'],
  # [0.05; 0.05 * view.CVaR95], [ones(1, J); X'], [1; view.mu]);
  p_ <- cbind(p_, tempEntropy$p_)
  # when maximizing, objective turns negative. Bringing it back to normal
  KLdiv <- rbind(KLdiv, -tempEntropy$optimizationPerformance$ml)

  # if change in KLdiv less than one percent, stop
  if(abs((KLdiv[i] - KLdiv[i - 1]) / KLdiv[i - 1])  < 0.01)
    doStop <- 1
}

# Display results

N <- length(s)
dev.new()
par(mfrow = c(2,1))
plot(1:N, KLdiv, type = "l", xlab = "s", ylab = "KL divergence")
lines(repmat(1:N, N, 1), repmat(t(KLdiv), N, 1), type = "p")
x <- seq(min(X), max(X), len = J)
tmp <- market.pdf(x)
tmp <- tmp / sum(tmp)
plot(t(X), tmp, type = "l", col = "blue", lwd = 1.5, xlab="",
     ylab = "", ylim = c(0, 0.006))
lines(t(X), p_[, ncol(p_)], type = "l", col = "red", lwd = 1.5)
lines(market.CVaR95, 0.0,  type = "p", pch = 17, col = "blue")
lines(view.CVaR95, 0.0,  type = "p", pch = 17, col = "red")
legend("topright", 1.9, c("Prior", "Posterior"),
       col = c("blue", "red"), lty = 1)

# zoom here
dev.new()
plot(t(X), tmp, "l", col = "blue", xlim = c(-4, -2), ylim = c(0, 0.001),
     xlab="", ylab = "")
lines(t(X), p_[, 1], col = "magenta")
lines(t(X), p_[, 2], col = "green")
lines(t(X), p_[, ncol(p_)], col = "red")
lines(market.CVaR95, 0.0,  type = "p", pch = 17, col = "blue")
lines(view.CVaR95, 0.0,  type = "p", pch = 17, col = "red")
