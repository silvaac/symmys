function CE = Satisfaction(Allocation,Market,InvestorProfile)

CE = Allocation'*diag(Market.CurrentPrices)*(1+Market.LinRets_EV) - ...
  1/(2*InvestorProfile.RiskPropensity)*Allocation'*diag(Market.CurrentPrices)*Market.LinRets_Cov*diag(Market.CurrentPrices)*Allocation ;