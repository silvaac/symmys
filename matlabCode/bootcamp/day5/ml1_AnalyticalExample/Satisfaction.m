function [CE,ExpectedValue,Variance] = Satisfaction(Allocation,Market,InvestorProfile)

ExpectedValues=diag(Market.CurrentPrices)*(1+Market.LinRets_EV);
Covariance=diag(Market.CurrentPrices)*Market.LinRets_Cov*diag(Market.CurrentPrices);

ExpectedValue=Allocation'*ExpectedValues;
Variance=Allocation'*Covariance*Allocation;
CE =  ExpectedValue - Variance/(2*InvestorProfile.RiskPropensity);