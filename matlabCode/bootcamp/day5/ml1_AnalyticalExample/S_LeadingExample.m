% this script performs the computations for the analytical example discussed
% in Sec. 6.1 of "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input investor's parameters
InvestorProfile.Budget=10000;
InvestorProfile.RiskPropensity=3000;
InvestorProfile.Confidence=.9;
InvestorProfile.BudgetAtRisk=.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input market parameters
NumAssets=5;
Min_Std=.06; Max_Std=.36; Step=(Max_Std-Min_Std)/(NumAssets-1);
Market.St_Devations  = [Min_Std : Step : Max_Std]';               
Market.LinRets_EV= .5*Market.St_Devations;                     
Market.Correlation  = .4;
Market.CurrentPrices=10*ones(NumAssets,1);                     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Correlation = (1-Market.Correlation) * eye(NumAssets) + Market.Correlation * ones(NumAssets,NumAssets);
Market.LinRets_Cov = diag(Market.St_Devations)*Correlation*diag(Market.St_Devations);

% compute optimal allocation, 
Optimal = ChoiceOptimal(Market,InvestorProfile);
[Optimal_CE,Optimal_ExpectedValue,Optimal_Variance] = Satisfaction(Optimal,Market,InvestorProfile);

% compute sub-optimal allocation
SubOptimal = ChoiceSubOptimal(Market,InvestorProfile);    
[SubOptimal_CE,SubOptimal_ExpectedValue,SubOptimal_Variance] = Satisfaction(SubOptimal,Market,InvestorProfile);

% compute frontier of feasible set
E=12000;  % expected value limit
[E,V,MV_ExpectedValue,MV_Variance,SR_ExpectedValue,SR_Variance]=ComputeFrontier(Market,InvestorProfile,E);

% compute var-constraint line
Var_ConstraintLine = (1-InvestorProfile.BudgetAtRisk)*InvestorProfile.Budget + sqrt(2*V)*erfinv(2*InvestorProfile.Confidence-1);

% compute iso-satisfaction lines
Optimal_SatisfactionLine = Optimal_CE + V/(2*InvestorProfile.RiskPropensity);
SubOptimal_SatisfactionLine = SubOptimal_CE + V/(2*InvestorProfile.RiskPropensity);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

% consraints
figure 

%frontier
Vertical=E(end:-1:1);
Horizontal=sqrt(V(end))+(max(sqrt(V))-min(sqrt(V)))/15*sin(2*pi*7*(Vertical-Vertical(end))/Vertical(end)-Vertical(1));
Fill_y=[E Vertical];
Fill_x=[sqrt(V) Horizontal];
h=fill(Fill_x,Fill_y,[.85 .85 .85])
set(h,'LineStyle','--')

Top_H=[min(sqrt(V)) : (max(sqrt(V))-min(sqrt(V)))/100 :  max(sqrt(V))];
Top=max(Var_ConstraintLine);
Top_V=Top+(max(E)-min(E))/15*sin(2*pi*3.5*(Top_H-min(sqrt(V)))/(max(sqrt(V))-min(sqrt(V)))  );
Left_V=[min(Var_ConstraintLine) : (max(Var_ConstraintLine)-min(Var_ConstraintLine))/100 : max(Var_ConstraintLine)];
Left_H=min(sqrt(V))+(max(sqrt(V))-min(sqrt(V)))/15*sin(2*pi*3.5*(Left_V-Left_V(end))/(Left_V(end)-Left_V(1)));


Fill_y=[Top_V max(Var_ConstraintLine)  min(Var_ConstraintLine) Left_V];
Fill_x=[Top_H max(sqrt(V)) min(sqrt(V)) Left_H];
hold on
h=fill(Fill_x,Fill_y,[.85 .85 .85])
set(h,'LineStyle','--')

[a,Index]=find(E>=Var_ConstraintLine);
Fill_y=E(Index);
Fill_x=sqrt(V(Index));
hold on
h=fill(Fill_x,Fill_y,[.7 .7 .7])

hold on
h=plot(sqrt(V),E);
set(h,'linewidth',1,'color','k');
hold on
h=plot(sqrt(V),Var_ConstraintLine);
set(h,'linewidth',1,'color','k');
hold on
h=plot(sqrt(V(Index)),E(Index));
set(h,'linewidth',2,'color','k');
hold on
Bottom=min(E(Index));Top=max(E(Index));
[a,Indexx]=find( (Var_ConstraintLine<=Top)  &  (Var_ConstraintLine>=Bottom) ); 
Index=intersect(Index,Indexx);
h=plot(sqrt(V(Index)),Var_ConstraintLine(Index));
set(h,'linewidth',2,'color','k');

grid on

% optimality
figure 

%frontier
Vertical=E(end:-1:1);
Horizontal=V(end)+(max(V)-min(V))/15*sin(2*pi*7*(Vertical-Vertical(1))/Vertical(end)-Vertical(1));
Fill_y=[E Vertical];
Fill_x=[V Horizontal];
h=fill(Fill_x,Fill_y,[.85 .85 .85]);
set(h,'LineStyle','--')
hold on
h=plot(V,E);
set(h,'linewidth',1,'color','k');
hold on
h=plot(V,Var_ConstraintLine,'LineStyle','--');
set(h,'linewidth',1,'color','k');

hold on  
h=plot(Optimal_Variance,Optimal_ExpectedValue,'.');
set(h,'markersize',20,'color','k');
hold on
h=plot(SubOptimal_Variance,SubOptimal_ExpectedValue,'.');
set(h,'markersize',20,'color','k');
hold on
[a,Index]=find(Optimal_SatisfactionLine<=max(E)); 
h=plot(V(Index),Optimal_SatisfactionLine(Index));
set(h,'linewidth',2,'color','k');
hold on
[a,Index]=find(SubOptimal_SatisfactionLine<=max(E)); 
h=plot(V(Index),SubOptimal_SatisfactionLine(Index));
set(h,'linewidth',2,'color','k');
grid on

