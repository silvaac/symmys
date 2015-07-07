% this script compares the EVT and sample estimates of the VaR with the true analytical VaR 
% under t-distribution assumptions
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci
% see also the technical appendices at symmys.com>book>downloads

clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs 

% market parameters (student t distribution)
mu=0;
sigma=1;
nu=100;

% EVT threshold
th=.95; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% confidence range for quantiles
c=[th :.001:.999];  

% generate simulations
NumSimul=100;
Psi =mu+sigma*trnd(nu,NumSimul,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical VaR
Q_an=mu+sigma*tinv(1-c,nu);

% empirical VaR
Q_em=prctile(Psi,(1-c)*100);

% EVT approximation
ThresholdIndex=ceil((1-th)*NumSimul); % percentile index
[Sort_Psi,Index] = sort(Psi);
psi_hat = Sort_Psi(ThresholdIndex);   % threshold level
Fpsi_hat = 1-th;                      % (empirical) cdf at threshold level 

Excess=psi_hat-Psi(Psi<psi_hat);
paramEsts=gpfit(Excess);
xi=paramEsts(1);
v=paramEsts(2);

Q_EVT = psi_hat + v/xi*(1-((1-c)/Fpsi_hat).^(-xi) );


figure
h1=plot(c,Q_an,'b');
hold on
h2=plot(c,Q_em,'g');
hold on
h3=plot(c,Q_EVT,'r');
legend([h1 h2 h3],'true','empirical','EVT',3)
grid on