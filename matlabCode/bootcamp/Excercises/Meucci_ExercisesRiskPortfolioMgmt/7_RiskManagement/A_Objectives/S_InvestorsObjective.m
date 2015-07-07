% This script familiarizes the users with the objectives of different investors
% in a highly non-normal bi-variate market of securities
% See Sec.5.1 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters of first marginal
nu_1=3;
s2_1=3;

% parameters of second marginal
mu_2=.1;
s2_2=.2;

% correlation in normal copula
r=.5;

% number of simulations
J=10000;  

% portfolio allocation
a=[1 2]';
% benchmark allocation
b=[2 1]'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute current prices
p_1=nu_1*s2_1;
p_2=exp(mu_2+.5*s2_2^2);
p=[p_1 p_2]';

% generate samnple of prices at the investment horizon
N=mvnrnd([0 0]',[1 r;r 1],J);
N_1=N(:,1);
N_2=N(:,2);

U_1=normcdf(N_1);
U_2=normcdf(N_2);

aa=nu_1/2;
bb=2*s2_1;
P_1=gaminv(U_1,aa,bb);
P_2=logninv(U_2,mu_2,sqrt(s2_2));

P=[P_1 P_2];

% generate sample of final wealth
W=P*a;

% generate sample of PnL
PnL=(P-ones(J,1)*p')*a;

% generate sample of benchmark-relative wealth
K=eye(2)-p*b'/(b'*p);
WRel=P*K'*a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
NumBins=round(10*log(J));
figure
plot(P_1,P_2,'.')
grid on
xlabel('P_1')
ylabel('P_2')

figure
hist(W,NumBins)
title('final wealth')
grid on

figure
hist(PnL,NumBins)
title('P&L')
grid on

figure
hist(WRel,NumBins)
title('benchmark-relative wealth')
grid on