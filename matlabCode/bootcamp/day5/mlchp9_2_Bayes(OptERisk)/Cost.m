function C_Plus = Cost(Allocation,Market,InvestorProfile)


aXi=Allocation'*diag(Market.CurrentPrices)*(1+Market.LinRets_EV);
aPhia=Allocation'*diag(Market.CurrentPrices)*Market.LinRets_Cov*diag(Market.CurrentPrices)*Allocation ;

C=(1-InvestorProfile.BaR)*InvestorProfile.Budget ...
        -aXi + sqrt(2*aPhia)*erfinv(2*InvestorProfile.Confidence-1);
C_Plus=max(C,0);    
