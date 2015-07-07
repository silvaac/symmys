% this script performs all the computations for a two-step optimal allocation, 
% where the first step is mean-variance and the second step is exact satisfaction computation
% see case study in Sec. 6.7 of "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Budget=10;
g=-9;                         % risk aversion parameter
Horizon=1;

Num_Scenarios=100000;
NumPortf=30;
Generate=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input data
if Generate   
  Num_Securities=8;
  Num_Observations=52*5;

  Rets_Correlation=eye(Num_Securities);
  Max_Ret=.12; Min_Ret=.03; Step=(Max_Ret-Min_Ret)/(Num_Securities-1);
  Rets_ExpectedValue=[Min_Ret : Step : Max_Ret]';
  Rets_StDeviation=2*Rets_ExpectedValue;
  Rets_Covariance=diag(Rets_StDeviation)*Rets_Correlation*diag(Rets_StDeviation);
  C=mvnrnd(1/52*Rets_ExpectedValue,1/52*Rets_Covariance,Num_Observations);
  P_0=rand(Num_Securities,1);
  Prices_TimeSeries=exp(cumsum(C));
  save Series_DB Prices_TimeSeries
else 
   load Series_DB;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate weekly invariants (compounded returns)
Week_Comp_Rets=log(Prices_TimeSeries(2:end,:))-log(Prices_TimeSeries(1:end-1,:));
Num_Observations=size(Week_Comp_Rets,1);
Num_Securities=size(Week_Comp_Rets,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% project market to horizon
Exp_Yr_Comp_Rets_Hat=mean(Week_Comp_Rets)'*52;
Cov_Yr_Comp_Rets_Hat=cov(Week_Comp_Rets)*52;
Current_Prices=Prices_TimeSeries(end,:)';
[Exp_Prices,Cov_Prices] = MktProjection(Exp_Yr_Comp_Rets_Hat,Cov_Yr_Comp_Rets_Hat,Current_Prices,Horizon);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% portolio optimization

% step 1: MV quadratic optimization to determine one-parameter frontier of quasi-optimal solutions ...
[ExpectedValue,Std_Deviation, Allocations] = ... % frontier in terms of number of shares
    EfficientFrontier(NumPortf, Cov_Prices, Exp_Prices,Current_Prices,Budget);
Rel_Allocations=Allocations.*(ones(NumPortf,1)*Current_Prices')/Budget; % frontier in terms of relative weights

% step 2: evaluate satisfaction for all allocations on the frontier ...
Returns_Scenarios=mvnrnd(Exp_Yr_Comp_Rets_Hat*Horizon,Cov_Yr_Comp_Rets_Hat*Horizon,Num_Scenarios);
Market_Scenarios=(ones(Num_Scenarios,1)*Current_Prices').*exp(Returns_Scenarios);
Store_Satisfaction=[];
for n=1:NumPortf
  Allocation=Allocations(n,:)';
  Objective_Scenario=Market_Scenarios*Allocation;
  Utility=(Objective_Scenario.^g)/g;
  Satisfaction=(g*mean(Utility))^(1/g);
  Store_Satisfaction=[Store_Satisfaction Satisfaction];
end

% ... and pick the best
[a,Optimal_Index]=max(Store_Satisfaction);
Optimal_Allocation=Allocations(Optimal_Index,:).*Current_Prices';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure

subplot(3,1,1) % MV frontier in moments coordinates
h=plot(Std_Deviation,ExpectedValue);
set(h,'linewidth',2,'color','k');
set(gca,'xlim',[Std_Deviation(1) Std_Deviation(end)],'ylim',[min(ExpectedValue) max(ExpectedValue)]);
grid on
xlabel('std. deviation')
ylabel('expected value')

subplot(3,1,2) % relative portfolio composition 
Data=cumsum(Rel_Allocations,2);
for n=1:Num_Securities
    x=[Std_Deviation(1); Std_Deviation; Std_Deviation(end)];
    y=[0; Data(:,Num_Securities-n+1); 0];
    hold on
    h=fill(x,y,[.85 .85 .85]-mod(n,2)*[.2 .2 .2]);
end
set(gca,'xlim',[Std_Deviation(1) Std_Deviation(end)],'ylim',[0 max(max(Data))])%,'Position',[0.13 0.1 0.775 0.45])
grid on
xlabel('std. deviation')
ylabel('portfolio composition')

subplot(3,1,3) % satisfaction as function of st.deviation on the frontier
h=plot(Std_Deviation,Store_Satisfaction);
set(h,'linewidth',2,'color','k');
set(gca,'xlim',[Std_Deviation(1) Std_Deviation(end)],'ylim',[min(Store_Satisfaction) max(Store_Satisfaction)]);
grid on
xlabel('std. deviation')
ylabel('satisfaction')