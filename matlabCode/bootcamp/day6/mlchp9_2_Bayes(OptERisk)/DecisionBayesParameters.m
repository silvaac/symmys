function Allocation = DecisionBayesParameters(Market,InvestorProfile,Prior)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate market parameters

Exp_LinRets_Hat=mean(Market.LinRetsSeries)';
Cov_LinRets_Hat=cov(Market.LinRetsSeries);
T=Market.T;

T_1=Prior.T_0+T;
Exp_LinRets_Post=(Prior.LinRets_EV*Prior.T_0+Exp_LinRets_Hat*T)/T_1;
Nu_1=Prior.Nu_0+T;
Cov_LinRets_Post=(Prior.LinRets_Cov*Prior.Nu_0+Cov_LinRets_Hat*T+(Exp_LinRets_Hat-Prior.LinRets_EV)*(Exp_LinRets_Hat-Prior.LinRets_EV)'*Prior.T_0*T/T_1)/Nu_1;

Exp_Prices_Hat=diag(Market.CurrentPrices)*(1+Exp_LinRets_Post);
Cov_Prices_Hat=diag(Market.CurrentPrices)*Cov_LinRets_Post*diag(Market.CurrentPrices);

% compute allocation
S=inv(Cov_Prices_Hat);
A_Hat=Market.CurrentPrices'*S*Market.CurrentPrices; 
B_Hat=Market.CurrentPrices'*S*Exp_Prices_Hat; 

Gamma = (InvestorProfile.Budget - InvestorProfile.RiskPropensity*B_Hat)/A_Hat;
Allocation = InvestorProfile.RiskPropensity*S*Exp_Prices_Hat + Gamma*S*Market.CurrentPrices;





