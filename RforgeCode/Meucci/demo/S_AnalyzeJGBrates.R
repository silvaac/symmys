#' This script applies the Inverse Call Transformation to Japanese government
#' rates and returns shadow rates, comparing them with rates and log-rates.
#'
#' @references
#' A. Meucci, A. Loregian - "Neither Normal not Lognormal: Modeling Interest
#' Rates Across all Regimes" \url{http://symmys.com/node/601}
#'
#' See Meucci's script "S_AnalyzeJGBrates.m"
#
#' @author Xavier Valls \email{xaviervallspla@@gmail.com}

data(JGB)

## process data

# time to maturity [1 2 3 5 7 10 15 20] years
tau <- JGB$TimeToMat[c( 1, 2, 3, 5, 7, 10, 11, 12)] 

# rates
y <- t(JGB$Yields[, c( 1, 2, 3, 5, 7, 10, 11, 12)]) 
Date <- JGB$date 
TimeStep <- 5  # daily (1) / weekly (5) observations
date <- Date[ seq(1,length(Date), TimeStep) ]
y<-y[ , seq(1,length(Date), TimeStep)]

# shadow rates, via Inverse Call Transformation
eta <-0.005
zeta <-0
icy <- InverseCallTransformation(y, tau, eta, zeta)

# log-rates
lny <- log(y)

# changes
dy <- diff(y)
dlny <- diff(log(y))
dx <- diff(icy)

t_ <- dim(dy)[2]
datetick <- seq(50,t_, 140)

## VISUAL
## rates/log-rates/shadow-rates plots
dev.new()
layout(matrix(c(1, 1, 2, 2, 3, 3, 4, 4), 4, 2, byrow = T), heights = c(3,3,3,3),
	   T)

matplot(t(y), type = "l", lty = 1, main = "Rates", xaxt = "n", ylab = "",
		col = topo.colors(20))
axis(1, at = datetick, labels = format(date[datetick], "%d-%b-%Y"))

matplot(t(lny), type = "l", lty = 1, main = "Log-rates", xaxt = "n", ylab = "",
		col = topo.colors(20)) 
axis(1, at = datetick, labels = format(date[datetick], "%d-%b-%Y"))

matplot(t(icy), type = "l", lty = 1, main = "Shadow rates", xaxt = "n", 
		ylab = "", col = topo.colors(20))
axis(1, at = datetick, labels = format(date[datetick], "%d-%b-%Y"))

plot.new()
par(mar=c(5.1,0,4.1, 0))
legend("bottom", y = c("1 year", "2 years", "3 years", "5 years", "7 years", 
	   "10 years", "15 years", "20 years"), lty = 1, ncol = 8,
	   col = topo.colors(20))
