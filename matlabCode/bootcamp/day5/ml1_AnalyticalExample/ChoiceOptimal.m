function Allocation = ChoiceOptimal(Market,InvestorProfile)


Exp_Prices=diag(Market.CurrentPrices)*(1+Market.LinRets_EV);
Cov_Prices=diag(Market.CurrentPrices)*Market.LinRets_Cov*diag(Market.CurrentPrices);

S=inv(Cov_Prices);
A=Market.CurrentPrices'*S*Market.CurrentPrices; 
B=Market.CurrentPrices'*S*Exp_Prices; 

Gamma = (InvestorProfile.Budget - InvestorProfile.RiskPropensity*B)/A;
Allocation = InvestorProfile.RiskPropensity*S*Exp_Prices + Gamma*S*Market.CurrentPrices;







