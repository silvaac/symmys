% this script search for cointegrated stat-arb strategies among swap contracts 
% see A. Meucci (2009) 
% "Review of Statistical Arbitrage, Cointegration, and Multivariate Ornstein-Uhlenbeck"
% available at ssrn.com

% Code by A. Meucci, April 2009
% Most recent version available at www.symmys.com > Teaching > MATLAB

clear; clc; close all
load DB_SwapParRates

%% estimation 
S=cov(Rates);
[E,Lam]=pcacov(S);

Y=Rates*E(:,end);

T=length(Y);
plot(Y)

coeff=regress(Y(2:end),[ones(T-1,1) Y(1:end-1)]);
a=coeff(1);
b=coeff(2);

theta=-log(b);
m=a/(1-b);

residuals=Y(2:end)-a-b*Y(1:end-1);
sig2=var(residuals)*2*theta/(1-exp(-2*theta));
sd=sqrt(sig2/(2*theta));

plot(Y)
hold on
plot(m*ones(1,T),'g')
hold on
plot(m*ones(1,T)+sd,'r')
hold on
plot(m*ones(1,T)-sd,'r')