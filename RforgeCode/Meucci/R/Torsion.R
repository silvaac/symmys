
#' @title Computes the Principal Components torsion and the Minimum Torsion
#'
#' @description  Computes the Principal Components torsion and the Minimum 
#' Torsion for diversification analysis, as described in A. Meucci, 
#' A. Santangelo, R. Deguest - "Measuring Portfolio Diversification Based on 
#' Optimized Uncorrelated Factors"
#'
#' @param  Sigma         [matrix] (n_ x n_) covariance matrix
#' @param  model         [string] choose between 'pca' and 'minimum-torsion' 
#'                                model
#' @param  method        [string] choose between 'approximate' and 'exact' 
#'                                method for 'minimum-torsion' model
#' @param  max_niter     [scalar] choose number of iterations of numerical 
#'                                algorithm 
#'  
#' @return t             [matrix] (n_ x n_) torsion matrix
#'
#' @references
#' A. Meucci, A. Santangelo, R. Deguest - "Measuring Portfolio Diversification 
#' Based on Optimized Uncorrelated Factors" \url{http://symmys.com/node/599}
#'
#' See Meucci's script "torsion.m"
#
#' @author Xavier Valls \email{xaviervallspla@@gmail.com}
#' @export

Torsion = function(Sigma, model, method = NULL, max_niter = 10000 ) {

  if (model == "pca") {
    # PCA decomposition
    e <- eigen(Sigma)
    lambda <- e$values
    ev <- e$vectors
    flip <- ev[1,] < 0
    # fix the sign of the eigenvector based on the sign of its first entry
    ev[, flip] <- -ev[, flip]
    index <- order(-lambda)

    # PCA torsion
    t <- t(ev[, index])

  } else if (model == "minimum-torsion") {
    # Correlation matrix
    sigma <- diag(Sigma) ^ (1 / 2)
    C <- diag(1 / sigma) %*% Sigma %*% diag(1 / sigma)
    c <- sqrtm(C)$B  # Riccati root of C

    if (method == "approximate") {
      t <- mrdivide(diag(sigma), c) %*% diag(1 / sigma)
    } else if (method == "exact") {
      n_ <- nrow(Sigma)

      # initialize
      d <- array(1, n_)
      f <- array(0, max_niter)
      for (i in 1:max_niter) {

        U <- diag(d) %*% c %*% c %*% diag(d)
        u <- sqrtm(U)$B
        q <- mldivide(u, (diag(d) %*% c))
        d <- diag(q %*% c)
        pi_ <- diag(d) %*% q # perturbation
        f[i] <- norm(c - pi_, type = "f")

        if(i > 1 && abs(f[i] - f[i - 1]) / f[i] / n_ <= 10 ^ (-8)) {
          f <- f[1:i]
          break
        } else if( i == max_niter && abs(f[max_niter] -
                  f[max_niter - 1]) / f[max_niter] / n_ > 10 ^ -8 ) {
          print(paste("number of max iterations reached: n_iter =", max_niter))
        }
      }
      x <- pi_ %*% solve(c)
      t <- diag(sigma) %*% x %*% diag(1/sigma)
    }
  }
  return(t)
}
