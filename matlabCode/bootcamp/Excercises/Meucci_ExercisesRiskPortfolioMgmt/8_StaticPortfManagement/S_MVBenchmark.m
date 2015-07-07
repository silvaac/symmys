% This script projects the distribution of the market invariants for the bond and stock markets 
% (i.e. the changes in yield to maturity and compounded returns) 
% from the estimation interval to the investment horizon 
% Then it computes the distribution of prices at the investment horizon and
% translates this distribution into the returns distribution
% Finally, it computes the mean-variance efficient frontier both for a
% total-return and for a benchmark-driven investor
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
tau=4/52;        % time to horizon expressed in years
tau_tilde=1/52;  % estimation period expressed in years

FlatCurve=.04;   
TimeToMat=5/52; % time to maturity of bond expressed in years

% parameters of the distribution of the changes in yield to maturity
u_minus_tau=TimeToMat-tau;
mu=0*u_minus_tau;
sigma=(20+5/4*u_minus_tau)/10000;

Num_Scenarios=100000;
Budget=100;
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
BondCurrent_Prices=exp(-FlatCurve*TimeToMat);

% project bond market to horizon

% generate changes in yield-to-maturity
DY_Scenarios = normrnd(mu*tau/tau_tilde,sigma*sqrt(tau/tau_tilde),Num_Scenarios,1); 
% compute the horizon prices, (3.81) in "Risk and Asset Allocation" - Springer
X=-u_minus_tau*DY_Scenarios;
BondMarket_Scenarios=BondCurrent_Prices_Shifted*exp(X); 

% MV inputs - analytical
Exp_Hrzn_DY_Hat=mu*tau/tau_tilde;
SDev_Hrzn_DY_Hat=sigma*sqrt(tau/tau_tilde);
Cov_Hrzn_DY_Hat=diag(SDev_Hrzn_DY_Hat)*diag(SDev_Hrzn_DY_Hat);
[BondExp_Prices,BondCov_Prices]=Dy2Prices(Exp_Hrzn_DY_Hat,Cov_Hrzn_DY_Hat,u_minus_tau,BondCurrent_Prices_Shifted)

% MV inputs - numerical
BondExp_Prices=mean(BondMarket_Scenarios)'
BondCov_Prices=cov(BondMarket_Scenarios)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put market together and compute returns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Current_Prices=[StockCurrent_Prices; BondCurrent_Prices];
Prices_Scenarios=[StockMarket_Scenarios BondMarket_Scenarios];
Rets_Scenarios=Prices_Scenarios./(ones(Num_Scenarios,1)*Current_Prices')-1;
E = mean(Rets_Scenarios)';
S = cov(Rets_Scenarios)

N=size(Prices_TimeSeries,2)+1;
w_b=ones(N,1)/N; % relative benchmar weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% portolio optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MV total return quadratic optimization to determine one-parameter frontier of quasi-optimal solutions 
NumPortf=40;
[ExpectedValue,Std_Deviation, Weights] = EfficientFrontierQPRets(NumPortf, S, E);
for k=1:NumPortf
    Rel_ExpectedValue(k)=(Weights(k,:)'-w_b)'*E;
    Rel_Std_Deviation(k)=sqrt((Weights(k,:)'-w_b)'*S*(Weights(k,:)'-w_b));
end

% MV benchmark-relative quadratic optimization to determine one-parameter frontier of quasi-optimal solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% benchmark-relative statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ExpectedValue_b,Std_Deviation_b, Weights_b] = EfficientFrontierQPRetsBench(NumPortf, S, E, w_b);
for k=1:NumPortf
    Rel_ExpectedValue_b(k)=(Weights_b(k,:)'-w_b)'*E;
    Rel_Std_Deviation_b(k)=sqrt((Weights_b(k,:)'-w_b)'*S*(Weights_b(k,:)'-w_b));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% frontiers in total return space
figure
h=plot(Std_Deviation,ExpectedValue);
set(h,'linewidth',2,'color','b');
hold on
h=plot(Std_Deviation_b,ExpectedValue_b);
set(h,'linewidth',2,'color','r');
legend('total ret','relative','location','best')
set(gca,'xlim',[Std_Deviation_b(1) Std_Deviation_b(end)],'ylim',[min(ExpectedValue_b) max(ExpectedValue_b)]);
xlabel('st.dev. rets.')
ylabel('exp.val rets.')

% frontiers in relative return space
figure
h=plot(Rel_Std_Deviation,Rel_ExpectedValue);
set(h,'linewidth',2,'color','b');
hold on
h=plot(Rel_Std_Deviation_b,Rel_ExpectedValue_b);
set(h,'linewidth',2,'color','r');
legend('total ret','relative','location','best')
set(gca,'xlim',[Rel_Std_Deviation_b(1) Rel_Std_Deviation_b(end)],'ylim',[min(Rel_ExpectedValue_b) max(Rel_ExpectedValue_b)]);
xlabel('TE rets.')
ylabel('EOP rets.')