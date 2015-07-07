function [Exp_Prices,Cov_Prices] = Dy2Prices(Exp_DY,Cov_DY,Times2Mat,CurrentPrices)
% this function computes the mean-variance inputs for zero-coupon bonds
% see (6.77)-(6.79) in "Risk and Asset Allocation"-Springer (2005), by A. Meucci


Mu=log(CurrentPrices)-Times2Mat.*Exp_DY;
Sigma=diag(Times2Mat.^2)*Cov_DY;

Exp_Prices=exp(Mu+(1/2)*diag(Sigma));
Cov_Prices=exp(Mu+(1/2)*diag(Sigma))*exp(Mu+(1/2)*diag(Sigma))'.*(exp(Sigma)-1);
