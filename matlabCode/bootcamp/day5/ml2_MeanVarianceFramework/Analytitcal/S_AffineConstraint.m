% this script shows the computation behind the analytical solution of the
% MV approach, namely when the only constraint is affine (e.g. budget constraint)
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input investor's parameters
InvestorProfile.Budget=.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input market parameters
NumAssets=5;
Min_EV=.03; Max_EV=.18; Step=(Max_EV-Min_EV)/(NumAssets-1);
Market.LinRets_EV = [Min_EV : Step : Max_EV]';               % hidden
Market.St_Devations = 2*Market.LinRets_EV;  % hidden
Market.Correlation  = .4;
Market.CurrentPrices=ones(NumAssets,1);     % not hidden

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Correlation = (1-Market.Correlation) * eye(NumAssets)   + Market.Correlation * ones(NumAssets,NumAssets);
Market.LinRets_Cov = diag(Market.St_Devations)*Correlation*diag(Market.St_Devations);

% compute frontier
[E,V,MV_ExpectedValue,MV_Variance,SR_ExpectedValue,SR_Variance]=ComputeFrontier(Market,InvestorProfile);

% line through origin, minimum variance, and maximum Sharpe ratio portfolios in variance vs. expected-value coordinates
Bottom=MV_Variance-1*(SR_Variance-MV_Variance); Top=SR_Variance; Step=(Top-Bottom)/100;
ThroughLine_V=[Bottom : Step : Top];
ThroughLine_E = MV_ExpectedValue + (SR_ExpectedValue - MV_ExpectedValue)/(SR_Variance - MV_Variance)*(ThroughLine_V-MV_Variance);

% line through origin and maximum Sharpe ratio portfolios in standard deviation vs. expected-value coordinates
Bottom=sqrt(MV_Variance)-2*(sqrt(SR_Variance)-sqrt(MV_Variance)); Top=sqrt(SR_Variance); Step=(Top-Bottom)/100;
ThroughLine_S=[Bottom : Step : Top];

%ThroughLine_S=[0 : sqrt(max(V))/100 : sqrt(max(V))];
ThroughLine_E2 =  SR_ExpectedValue/sqrt(SR_Variance)*ThroughLine_S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures

figure % mean - variance

Vertical=E(end:-1:1); % frontier
Horizontal=V(end)+(max(V)-min(V))/15*sin(2*pi*5*(Vertical-Vertical(1))/(max(Vertical)-min(Vertical)));
Fill_y=[E Vertical];
Fill_x=[V Horizontal];
h=fill(Fill_x,Fill_y,[.85 .85 .85]);
set(h,'linestyle','--')
hold on
h=plot(V(find(E>=MV_ExpectedValue)),E(find(E>=MV_ExpectedValue)));
set(h,'linewidth',2,'color','k');

hold on
h=plot(MV_Variance,MV_ExpectedValue,'.');
set(h,'markersize',20,'color','k');
hold on
h=plot(SR_Variance,SR_ExpectedValue,'.');
set(h,'markersize',20,'color','k');
hold on
h=plot(ThroughLine_V,ThroughLine_E);
set(h,'linewidth',1,'color','k');

grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure % mean - standard deviation

Vertical=E(end:-1:1); % frontier
Horizontal=sqrt(V(end))+(max(sqrt(V))-min(sqrt(V)))/15*sin(2*pi*5*(Vertical-Vertical(1))/(max(Vertical)-min(Vertical)));
Fill_y=[E Vertical];
Fill_x=[sqrt(V) Horizontal];
h=fill(Fill_x,Fill_y,[.85 .85 .85]);
set(h,'linestyle','--');
hold on
h=plot(sqrt(V(find(E>=MV_ExpectedValue))),E(find(E>=MV_ExpectedValue)));
set(h,'linewidth',2,'color','k');

h=plot(sqrt(MV_Variance),MV_ExpectedValue,'.');
set(h,'markersize',20,'color','k');
hold on
h=plot(sqrt(SR_Variance),SR_ExpectedValue,'.');
set(h,'markersize',20,'color','k');
hold on
h=plot(ThroughLine_S,ThroughLine_E2);
set(h,'linewidth',1,'color','k');

grid on

