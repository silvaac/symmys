function [Exp_Prices,Cov_Prices] = MktProjection(Exp_Yr_Comp_Rets,Cov_Yr_Comp_Rets,Starting_Prices,Horizon)

Mu=Exp_Yr_Comp_Rets*Horizon;
Sigma=Cov_Yr_Comp_Rets*Horizon;

M=exp(Mu+(1/2)*diag(Sigma));
Exp_Prices=diag(Starting_Prices)*M;

S=exp(Mu+(1/2)*diag(Sigma))*exp(Mu+(1/2)*diag(Sigma))'.*(exp(Sigma)-1);
Cov_Prices=diag(Starting_Prices)*S*diag(Starting_Prices);