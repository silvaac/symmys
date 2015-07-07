% This script projects the distribution of the market invariants for the bond and stock markets 
% (i.e. the changes in yield to maturity and compounded returns) 
% from the estimation interval to the investment horizon 
% Then it computes the distribution of prices at the investment horizon 
% and performs the two-step mean-variance optimization 
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
tau=4/52;        % time to horizon expressed in years
tau_tilde=1/52;  % estimation period expressed in years

FlatCurve=.04;   
TimesToMat=[4 5 10 52 520]'/52; % time to maturity of selected bonds expressed in years

% parameters of the distribution of the changes in yield to maturity
u_minus_tau=TimesToMat-tau;
nu=8;
mus=0*u_minus_tau;
sigmas=(20+5/4*u_minus_tau)/10000;

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
StockExp_Prices=mean(StockMarket_Scenarios)'
StockCov_Prices=cov(StockMarket_Scenarios)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bond market projection to horizon and pricing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BondCurrent_Prices_Shifted=exp(-FlatCurve*u_minus_tau);
BondCurrent_Prices=exp(-FlatCurve*TimesToMat);

% project bond market to horizon
N=length(TimesToMat); % number of bonds

% generate common source of randomness
U=rand(Num_Scenarios,1);
BondMarket_Scenarios=zeros(Num_Scenarios,N);
for n=1:N
    % generate co-dependent changes in yield-to-maturity
    DY_Scenarios = norminv(U,mus(n)*tau/tau_tilde,sigmas(n)*sqrt(tau/tau_tilde)); 

    % compute the horizon prices, (3.81) in "Risk and Asset Allocation" - Springer
    X=-u_minus_tau(n)*DY_Scenarios;
    BondMarket_Scenarios(:,n)=BondCurrent_Prices_Shifted(n)*exp(X); 
end

% MV inputs - analytical
Exp_Hrzn_DY_Hat=mus*tau/tau_tilde;
SDev_Hrzn_DY_Hat=sigmas*sqrt(tau/tau_tilde);
Corr_Hrzn_DY_Hat=ones(N); % full co-dependence
Cov_Hrzn_DY_Hat=diag(SDev_Hrzn_DY_Hat)*Corr_Hrzn_DY_Hat*diag(SDev_Hrzn_DY_Hat);
[BondExp_Prices,BondCov_Prices]=Dy2Prices(Exp_Hrzn_DY_Hat,Cov_Hrzn_DY_Hat,u_minus_tau,BondCurrent_Prices_Shifted)

% MV inputs - numerical
BondExp_Prices=mean(BondMarket_Scenarios)'
BondCov_Prices=cov(BondMarket_Scenarios)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% portolio optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% step 1: MV quadratic optimization to determine one-parameter frontier of quasi-optimal solutions ...
E=[StockExp_Prices
    BondExp_Prices(2)];
S = blkdiag(StockCov_Prices,BondCov_Prices(2,2));
Current_Prices=[StockCurrent_Prices
    BondCurrent_Prices(2)];
Market_Scenarios=[StockMarket_Scenarios BondMarket_Scenarios(:,2)];

NumPortf=40;
[ExpectedValue,Std_Deviation, Allocations] = ... % frontier with QP (no short-sales)
    EfficientFrontierQPPrices(NumPortf, S, E,Current_Prices,Budget);

% step 2: ...evaluate satisfaction for all allocations on the frontier ...
Store_Satisfaction=[];
for n=1:NumPortf
  Allocation=Allocations(n,:)';
  Objective_Scenario=Market_Scenarios*Allocation;
  Utility=-exp(-1/Zeta*Objective_Scenario);
  ExpU=mean(Utility);  
  Satisfaction=-Zeta*log(-ExpU);
  Store_Satisfaction=[Store_Satisfaction Satisfaction];
end

% ... and pick the best
[a,Optimal_Index]=max(Store_Satisfaction);
Optimal_Allocation=Allocations(Optimal_Index,:)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

subplot(2,1,1) % prices MV frontier 
h=plot(Std_Deviation,ExpectedValue);
set(h,'linewidth',2,'color','k');
set(gca,'xlim',[Std_Deviation(1) Std_Deviation(end)],'ylim',[min(ExpectedValue) max(ExpectedValue)]);
xlabel('st.dev. prices')
ylabel('exp.val. prices')

subplot(2,1,2) % satisfaction as function of st.deviation on the frontier
h=plot(Std_Deviation,Store_Satisfaction);
set(h,'linewidth',2,'color','k');
set(gca,'xlim',[Std_Deviation(1) Std_Deviation(end)],'ylim',[min(Store_Satisfaction) max(Store_Satisfaction)]);
xlabel('st.dev. prices')
ylabel('satisfaction')
