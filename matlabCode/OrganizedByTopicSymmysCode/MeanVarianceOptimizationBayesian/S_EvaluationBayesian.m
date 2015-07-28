% this script evaluates the Bayesian allocation, which replaces the true, unknown market parameter 
% in the optimal allocation policy with a Bayesian classical-equivalent (point) estimate. 
% the script displays the distribution of satisfaction, cost of constraint violation and opportunity cost 
% for each value of the market stress-test parameters (in this case the correlation)
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment inputs
NumScenarios=5000; 
NumCorrelations=5; 
Bottom_Correlation=0;
Top_Correlation=.99;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input investor's parameters
InvestorProfile.Budget=10000;
InvestorProfile.RiskPropensity=3;
InvestorProfile.Confidence=.99;
InvestorProfile.BaR=.1; %.2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input market parameters
NumAssets=5;
a=.5;    % effect of correlation on expected values and volatility (hidden)
Bottom=.06; Top=.36; Step=(Top-Bottom)/(NumAssets-1); v=[Bottom : Step : Top]';  % volatility vector
Market.T=30;                                % not hidden
Market.CurrentPrices=10*ones(NumAssets,1);  % not hidden

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input prior in Bayes approach;
Rho=.5;        % market overall correlations
Prior.St_Devations = (1+a*Rho)*v;  
Prior.LinRets_EV = .5*Prior.St_Devations;       
Correlation = (1-Rho) * eye(NumAssets) + Rho * ones(NumAssets,NumAssets);
Prior.LinRets_Cov = diag(Prior.St_Devations)*Correlation*diag(Prior.St_Devations);
Prior.T_0=Market.T;
Prior.Nu_0=Market.T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumCases=5000;                % compute the empirical pdf of the prior correlation for the plot
av_cor = InvWishCorr(Prior.LinRets_Cov,Prior.Nu_0,NumCases); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Step=(Top_Correlation-Bottom_Correlation)/(NumCorrelations-1);
Overall_Correlations=[Bottom_Correlation : Step : Top_Correlation]; 

Suboptimal.Store_Satisfaction=[]; 
Suboptimal.Store_CostConstraints=[]; 
Suboptimal.Store_OppCost=[];
Optimal.Store_Satisfaction=[];
for t=1:length(Overall_Correlations)
 
  Cycles_to_go=length(Overall_Correlations)-t+1 % display some info on the main window screen to know what's going on

  % input the (hidden) market parameters (only correlations, we assume standard deviations and expected values fixed and known)
  Market.St_Devations = (1+a*Overall_Correlations(t))*v;  % hidden
  Market.LinRets_EV = .5*Market.St_Devations;       % hidden
  Correlation = (1-Overall_Correlations(t)) * eye(NumAssets) + Overall_Correlations(t) * ones(NumAssets,NumAssets);
  Market.LinRets_Cov = diag(Market.St_Devations)*Correlation*diag(Market.St_Devations);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute optimal allocation, only possible if hidden parameters were known: thus it is not a "decision", we call it a "choice"
  AllocationOpt = ChoiceOptimal(Market,InvestorProfile);
  Satisfaction_Optimal = Satisfaction(AllocationOpt,Market,InvestorProfile);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % choose allocation based on available information
    Store_TrueSatisfaction=[]; Store_CostConstraints=[]; 
  for s=1:NumScenarios
    
    Market.LinRetsSeries=mvnrnd(Market.LinRets_EV,Market.LinRets_Cov,Market.T); % generate one scenario i_T of information I_T
    
    Allocation = DecisionBayesParameters(Market,InvestorProfile,Prior);       % scenario-dependent decision: tries to pick the optimal 
    TrueSatisfaction = Satisfaction(Allocation,Market,InvestorProfile);       % solution by Bayes estimation of the paramters
    CostConstraints=Cost(Allocation,Market,InvestorProfile);                  % the estimation error is there, but tapered    

    Store_TrueSatisfaction = [Store_TrueSatisfaction TrueSatisfaction];
    Store_CostConstraints = [Store_CostConstraints CostConstraints];
    
  end
  
  Suboptimal.Store_CostConstraints=[Suboptimal.Store_CostConstraints; Store_CostConstraints];
  Suboptimal.Store_Satisfaction=[Suboptimal.Store_Satisfaction; Store_TrueSatisfaction];
  Suboptimal.Store_OppCost=[Suboptimal.Store_OppCost; Satisfaction_Optimal-Store_TrueSatisfaction+Store_CostConstraints];
  Optimal.Store_Satisfaction=[Optimal.Store_Satisfaction; Satisfaction_Optimal];
    
end

PlotEvaluationBayesian