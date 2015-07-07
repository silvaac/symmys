function [E,V,MV_ExpectedValue,MV_Variance,SR_ExpectedValue,SR_Variance]=ComputeFrontier(Market,InvestorProfile,E)

% compute useful parameters
ExpectedValues=diag(Market.CurrentPrices)*(1+Market.LinRets_EV);
Covariance=diag(Market.CurrentPrices)*Market.LinRets_Cov*diag(Market.CurrentPrices);
S=inv(Covariance);
A=Market.CurrentPrices'*S*Market.CurrentPrices; 
B=Market.CurrentPrices'*S*ExpectedValues; 

% compute minimum-variance and maximum Sharpe ratio portfolios and their coordinates
MV_Portf=InvestorProfile.Budget*S*Market.CurrentPrices/A;
MV_ExpectedValue=MV_Portf'*ExpectedValues; 
MV_Variance=MV_Portf'*Covariance*MV_Portf; 
SR_Portf=InvestorProfile.Budget*S*ExpectedValues/B;
SR_ExpectedValue=SR_Portf'*ExpectedValues; 
SR_Variance=SR_Portf'*Covariance*SR_Portf; 

% compute upper branch of the feasible set (i.e., the frontier)
Step=InvestorProfile.Budget/100;
Step=(E-MV_ExpectedValue)/100;
E=[2*MV_ExpectedValue-E : Step : E];
V=[];
for i=1:length(E);
  Curve_Portf=MV_Portf + (E(i)-MV_ExpectedValue)* (SR_Portf-MV_Portf)/(SR_ExpectedValue-MV_ExpectedValue);
  Curve_Variance=Curve_Portf'*Covariance*Curve_Portf;
  V=[V Curve_Variance];
end
