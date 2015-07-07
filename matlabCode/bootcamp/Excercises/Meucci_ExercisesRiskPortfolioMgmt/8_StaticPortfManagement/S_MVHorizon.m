% This script projects the distribution of the market invariants for stock market (i.e. compounded returns) 
% from the estimation interval to the investment horizon 
% Then it computes the distribution of prices at the investment horizon 
% and performs the two-step mean-variance optimization in terms of returns and relative portfolio weights
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
tau=1/252;        % time to horizon expressed in years
tau_tilde=1/52;  % estimation period expressed in years
Num_Scenarios=100000;
Budget=100;
Zeta=10; % risk aversion parameter
load StockSeries;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimation of weekly invariants stock market (compounded returns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Week_C=log(Prices_TimeSeries(2:end,:))-log(Prices_TimeSeries(1:end-1,:));
[T,N]=size(Week_C);

la=.1;
Shrk_Exp=zeros(N,1);
Exp_C_Hat=(1-la)*mean(Week_C)'+la*Shrk_Exp;

lb=.1;
Shrk_Cov=eye(N)*trace(cov(Week_C))/N;
Cov_C_Hat=(1-lb)*cov(Week_C)+lb*(Shrk_Cov);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stock market projection to horizon and pricing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Exp_Hrzn_C_Hat=Exp_C_Hat*tau/tau_tilde;
Cov_Hrzn_C_Hat=Cov_C_Hat*tau/tau_tilde;
StockCompReturns_Scenarios=mvnrnd(Exp_Hrzn_C_Hat,Cov_Hrzn_C_Hat,Num_Scenarios);

StockCurrent_Prices=Prices_TimeSeries(end,:)';
StockMarket_Scenarios=(ones(Num_Scenarios,1)*StockCurrent_Prices').*exp(StockCompReturns_Scenarios);

% MV inputs - analytical
[StockExp_Prices,StockCov_Prices] = ComRets2Prices(Exp_Hrzn_C_Hat,Cov_Hrzn_C_Hat,StockCurrent_Prices)
% MV inputs - numerical
StockExp_Prices=mean(StockMarket_Scenarios)';
StockCov_Prices=cov(StockMarket_Scenarios);

StockExp_LinRets=StockExp_Prices./StockCurrent_Prices-1;
StockCov_LinRets=diag(1./StockCurrent_Prices)*StockCov_Prices*diag(1./StockCurrent_Prices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% portolio optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% step 1: MV quadratic optimization to determine one-parameter frontier of quasi-optimal solutions ...
NumPortf=40;
[ExpectedValue,Std_Deviation, Weights] = EfficientFrontierQPRets(NumPortf, StockCov_LinRets, StockExp_LinRets);

% step 2: ...evaluate satisfaction for all allocations on the frontier ...
Store_Satisfaction=[];
for n=1:NumPortf
  Allocation=Weights(n,:)'*Budget./StockCurrent_Prices;
  Objective_Scenario=StockMarket_Scenarios*Allocation;
  Utility=-exp(-1/Zeta*Objective_Scenario);
  ExpU=mean(Utility);  
  Satisfaction=-Zeta*log(-ExpU);
  Store_Satisfaction=[Store_Satisfaction Satisfaction];
end

% ... and pick the best
[a,Optimal_Index]=max(Store_Satisfaction);
Optimal_Allocation=Weights(Optimal_Index,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

subplot(2,1,1) % rets MV frontier 
h=plot(Std_Deviation,ExpectedValue);
set(h,'linewidth',2,'color','k');
set(gca,'xlim',[Std_Deviation(1) Std_Deviation(end)],'ylim',[min(ExpectedValue) max(ExpectedValue)]);
xlabel('st.dev. rets.')
ylabel('exp.val rets.')

subplot(2,1,2) % satisfaction as function of st.deviation on the frontier
h=plot(Std_Deviation,Store_Satisfaction);
set(h,'linewidth',2,'color','k');
set(gca,'xlim',[Std_Deviation(1) Std_Deviation(end)],'ylim',[min(Store_Satisfaction) max(Store_Satisfaction)]);
xlabel('st.dev. rets')
ylabel('satisfaction')
