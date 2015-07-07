function [Expected_Value,Covariance,Standard_Deviation,Correlation]=LogNormalParam2Statistics(Mu,Sigma)
% this function computes analytically the summary statistics of a lognormal variable
% given the input parameters
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci



Expected_Value=exp(Mu+(1/2)*diag(Sigma));

Covariance=exp(Mu+(1/2)*diag(Sigma))*exp(Mu+(1/2)*diag(Sigma))'.*(exp(Sigma)-1);

Standard_Deviation=sqrt(diag(Covariance));

Correlation=diag(1./Standard_Deviation)*Covariance*diag(1./Standard_Deviation);
