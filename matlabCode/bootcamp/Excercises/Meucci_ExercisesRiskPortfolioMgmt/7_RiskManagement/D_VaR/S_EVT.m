% this script computes the quantile (VaR) 
% - analytically, under the Student t assumption for the objective
% - in simulations, using the sample quantile
% - using the extreme value theory approximation
% See Sec 5.5 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs 

% market parameters (student t distribution)
m=1;
s=2;
nu=7;

th=.95; % EVT threshold
c=[th : .001 : .999];  % confidence range for quantiles

NumSimul=200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical
Q_an=m+s*tinv(1-c,nu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulations

% generate objective's scenarios
X=trnd(nu,NumSimul/2,1);
X=[X         % symmetrize simulations
    -X];
Psi = m + s*X;
Q_simul=prctile(Psi,(1-c)*100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVT approximation
psi_hat = prctile(Psi,(1-th)*100);
Excess=psi_hat-Psi(Psi<psi_hat);
xi_v=gpfit(Excess);
xi=xi_v(1);
v=xi_v(2);

Fpsi_hat = 1-th;
Q_EVT = psi_hat + v/xi*(1-((1-c)/Fpsi_hat).^(-xi) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(c,Q_an)
hold on 
plot(c,Q_simul,'g')
hold on
plot(c,Q_EVT,'r')
legend('exact','simulations','EVT',3)
grid on
xlabel('Confidence, c')
ylabel('Quantile based satisfaction, Q_c(\alpha)')
