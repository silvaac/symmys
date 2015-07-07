% this script evaluates a generic allocation decision (in this case the 
% "best performer" strategy, which fully invest the budget in the last period's best performer)
% it displays the distribution of satisfaction, cost of constraint violation and opportunity cost 
% for each value of the market stress-test parameters (in this case the correlation)
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment inputs
NumScenarios=10000;  
NumCorrelations=5;
Bottom_Correlation=0;
Top_Correlation=.99;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input investor's parameters
InvestorProfile.Budget=10000;
InvestorProfile.RiskPropensity=30;
InvestorProfile.Confidence=.9;
InvestorProfile.BaR=.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input market parameters
NumAssets=10;  
a=.5;                                      % effect of correlation on expected values and volatility (hidden)
Bottom=.06; Top=.36; Step=(Top-Bottom)/(NumAssets-1); v=[Bottom : Step : Top]';                 % volatility vector
Market.T=20;                                % not hidden
Market.CurrentPrices=10*ones(NumAssets,1);  % not hidden

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Step=(Top_Correlation-Bottom_Correlation)/(NumCorrelations-1);
Overall_Correlations=[Bottom_Correlation : Step : Top_Correlation]; 

Suboptimal.StrsTst_Satisfaction=[]; 
Suboptimal.StrsTst_CostConstraints=[]; 
Suboptimal.StrsTst_OppCost=[];
Optimal.StrsTst_Satisfaction=[];
for t=1:length(Overall_Correlations)
 
  Cycles_to_go=length(Overall_Correlations)-t+1 % display some info on the main window screen to know what's going on

  % input the (hidden) market parameters (only correlations, we assume standard deviations and expected values fixed and known)
  Market.St_Devations = (1+a*Overall_Correlations(t))*v;  % hidden
  Market.LinRets_EV = .5*Market.St_Devations;       % hidden

  Correlation = (1-Overall_Correlations(t)) * eye(NumAssets) + Overall_Correlations(t) * ones(NumAssets,NumAssets);
  Market.LinRets_Cov = diag(Market.St_Devations)*Correlation*diag(Market.St_Devations);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute optimal allocation, only possible if hidden parameters were known: thus it is not a "decision", we call it a "choice"
  Allocation = ChoiceOptimal(Market,InvestorProfile);
  Satisfaction_Optimal = Satisfaction(Allocation,Market,InvestorProfile);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % choose allocation based on available information
  StrsTst_TrueSatisfaction=[]; StrsTst_CostConstraints=[];
  for s=1:NumScenarios
    
    Market.LinRetsSeries=mvnrnd(Market.LinRets_EV,Market.LinRets_Cov,Market.T); % generate scenarios i_T of information I_T
    
    Allocation = DecisionBestPerformer(Market,InvestorProfile);   % scenario-dependent decision that tries to pick the optimal allocation
    TrueSatisfaction = Satisfaction(Allocation,Market,InvestorProfile); 
    CostConstraints=Cost(Allocation,Market,InvestorProfile); 
    
    StrsTst_TrueSatisfaction = [StrsTst_TrueSatisfaction TrueSatisfaction];
    StrsTst_CostConstraints = [StrsTst_CostConstraints CostConstraints];
    
  end
  
  Suboptimal.StrsTst_CostConstraints=[Suboptimal.StrsTst_CostConstraints; StrsTst_CostConstraints];
  Suboptimal.StrsTst_Satisfaction=[Suboptimal.StrsTst_Satisfaction; StrsTst_TrueSatisfaction];
  Suboptimal.StrsTst_OppCost=[Suboptimal.StrsTst_OppCost; Satisfaction_Optimal-StrsTst_TrueSatisfaction+StrsTst_CostConstraints];
  Optimal.StrsTst_Satisfaction=[Optimal.StrsTst_Satisfaction; Satisfaction_Optimal];
    
end

PlotEvaluationGeneric