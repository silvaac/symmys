function [Exp_Prices,Cov_Prices] = ComRets2Prices(Exp_Comp_Rets,Cov_Comp_Rets,Starting_Prices)
% this function computes the mean-variance inputs for stocks
% see (6.77)-(6.79) in "Risk and Asset Allocation"-Springer (2005), by A. Meucci


Mu=log(Starting_Prices)+Exp_Comp_Rets;
Sigma=Cov_Comp_Rets;

Exp_Prices=exp(Mu+(1/2)*diag(Sigma));
Cov_Prices=exp(Mu+(1/2)*diag(Sigma))*exp(Mu+(1/2)*diag(Sigma))'.*(exp(Sigma)-1);