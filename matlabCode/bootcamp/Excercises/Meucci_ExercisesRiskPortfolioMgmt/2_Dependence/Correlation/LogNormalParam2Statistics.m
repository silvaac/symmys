function [Expected_Value,Covariance,Standard_Deviation,Correlation]=LogNormalParam2Statistics(Mu,Sigma)

Expected_Value=exp(Mu+(1/2)*diag(Sigma));

Covariance=exp(Mu+(1/2)*diag(Sigma))*exp(Mu+(1/2)*diag(Sigma))'.*(exp(Sigma)-1);

Standard_Deviation=sqrt(diag(Covariance));

Correlation=diag(1./Standard_Deviation)*Covariance*diag(1./Standard_Deviation);
