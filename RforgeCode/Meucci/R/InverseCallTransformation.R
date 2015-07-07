#' @title Computes the Inverse Call Transformation
#'
#' @description  Computes the Inverse Call Transformation and returns shadow 
#' rates, as described in A. Meucci, A. Loregian -"Neither Normal not Lognormal:
#' Modeling Interest Rates Across all Regimes".
#'
#' @param  Rates           [matrix] (N x length(timeseries)) number of portfolio
#' in the efficient frontier
#' @param  Tau             [vector] (N x 1) vector containing the times to
#' maturity corresponding to the rows of the rates matrix
#' @param  Eta             [scalar] Inverse-call transformation parameter
#' @param  Zeta            [scalar] Inverse-call transformation parameter
#'  
#' @return x               [vector] (NumPortf x 1) shadow rates, computed from
#' rates via inverse-call transformation
#'
#' @note                   Smoothing parameter s obtained as s=eta*exp(zeta*tau)
#'  
#' @references
#' A. Meucci, A. Loregian - "Neither Normal not Lognormal: Modeling Interest
#' Rates Across all Regimes" \url{http://symmys.com/node/601}
#'
#' See Meucci's script "InverseCallTransformation.m".
#
#' @author Xavier Valls \email{xaviervallspla@@gmail.com}
#' @export

InverseCallTransformation <- function(rates, tau, eta, zeta){

  t_ <- dim(rates)[2]
  x  <- matrix(0, nrow(rates), ncol(rates))
  s  <- eta * exp( zeta * tau )


  for (v in 1:length(tau)){
    # inverse call transformation

    x0 <- 0  # initialization

    CallFit <- function(tmpX, y = rates[v, t], sigma = s[v]){
      c <- BachelierCallPrice(tmpX, sigma)
      return(F = y - c)
    }
    for (t in 1:t_){
      # fitting inverse call
      x[v,t] <- lsqnonlin(fun = CallFit, x0 = x0)$x
    }
  }

  return(x)
}

# Bachelier call pricing function
# Call function (zero-strike call option price profile according to the
# Bachelier pricing function)
#
# s: smoothing parameter

BachelierCallPrice <- function( x , s){
  return(x * pnorm( x / s ) + s %*% dnorm(x / s))
}
