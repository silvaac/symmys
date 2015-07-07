% this script shows how overall correlation affects the maximum level of satisfaction, as measured 
% by the certainty equivalent of an exponential utility function in a multivariate normal market
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input investor's parameters
InvestorProfile.Budget=10000;
InvestorProfile.RiskPropensity=20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input market parameters
NumAssets=4;
Overall_Correlations=[0 :.01 : .99]; 
Min_EV=.03; Max_EV=.18; Step=(Max_EV-Min_EV)/(NumAssets-1);
Market.LinRets_EV = [Min_EV : Step : Max_EV]';       % hidden
Market.St_Devations = 2*Market.LinRets_EV;           % hidden
Market.CurrentPrices=10*ones(NumAssets,1);           % not hidden

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Store_CE=[]; 
for t=1:length(Overall_Correlations)
    
    % input the remaining market parameters 
    Correlation = (1-Overall_Correlations(t)) * eye(NumAssets) + Overall_Correlations(t) * ones(NumAssets,NumAssets);
    Market.LinRets_Cov = diag(Market.St_Devations)*Correlation*diag(Market.St_Devations);
    
    % compute optimal allocation
    Market.Exp_Prices=diag(Market.CurrentPrices)*(1+Market.LinRets_EV);
    Market.Cov_Prices=diag(Market.CurrentPrices)*Market.LinRets_Cov*diag(Market.CurrentPrices);
    
    S=inv(Market.Cov_Prices);
    A=Market.CurrentPrices'*S*Market.CurrentPrices; 
    B=Market.CurrentPrices'*S*Market.Exp_Prices; 
    Gamma = (InvestorProfile.Budget - InvestorProfile.RiskPropensity*B)/A;
    
    Allocation = InvestorProfile.RiskPropensity*S*Market.Exp_Prices + Gamma*S*Market.CurrentPrices;
    CE = Allocation'*Market.Exp_Prices -1/(2*InvestorProfile.RiskPropensity)*Allocation'*Market.Cov_Prices*Allocation;
       
    Store_CE=[Store_CE; CE];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure % optimal allocation vs. fixed allocation
h=plot(Overall_Correlations,Store_CE);
grid on
box off