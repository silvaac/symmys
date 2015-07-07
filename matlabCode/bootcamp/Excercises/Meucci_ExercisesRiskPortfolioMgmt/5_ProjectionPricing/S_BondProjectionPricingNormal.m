% This script projects the distribution of the market invariants for the bond markets 
% (i.e. the changes in yield to maturity) 
% from the estimation interval to the investment horizon 
% Then it computes the distribution of prices at the investment horizon 
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
tau=1/52;        % time to horizon expressed in years
tau_tilde=1/52;  % estimation period expressed in years

FlatCurve=.04;   
TimesToMat=[4 5 10 52 520]'/52; % time to maturity of selected bonds expressed in years

% parameters of the distribution of the changes in yield to maturity
u_minus_tau=TimesToMat-tau;
mus=0*u_minus_tau;
sigmas=(20+5/4*u_minus_tau)/10000;

Num_Scenarios=100000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bond market projection to horizon and pricing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BondCurrent_Prices_Shifted=exp(-FlatCurve*u_minus_tau);
BondCurrent_Prices=exp(-FlatCurve*TimesToMat);

% project bond market to horizon
N=length(TimesToMat); % number of bonds
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

