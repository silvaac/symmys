% this script computes the expected shortfall and the contributions to ES from each security 
% - analytically, under the elliptical-uniform assumption for the market
% - in simulations, using the conditional expectation definition of the contributions
% See Sec 5.5 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs 

% number of assets
N=10; 

% market parameters (uniform on ellipsoid)
Mu=rand(N,1);
A=rand(N,N)-.5;
Sigma=A*A';

% allocation
a=rand(N,1)-.5;

% quantile level
c=.95;

NumSimul=10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate market scenarios
X=GenerateUniform2(NumSimul,N); % uniform on sphere
M = ones(NumSimul,1)*Mu' + X*A';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulations

% compute and sort the objective
Psi = M*a;
Q_sim=prctile(Psi,(1-c)*100)

e=mean(abs(a))/100; % perturbation
for n=1:N
    % compute gradient

    a_e=a;
    a_e(n)=a(n)+e;
    
    Psi_e = M*a_e;
    Q_sim_e=prctile(Psi_e,(1-c)*100);
    DQ_simul(n)=(Q_sim_e-Q_sim)/e;

end
% compute contributions
ContrQ_simul=a.*DQ_simul';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical

% compute quantile of standardized marginal (1-dim generator) in simulations
% this does NOT depend on the allocation...
gc=prctile(X(:,1),(1-c)*100);

% ...the dependence on the allocation is analytical
Q_an=Mu'*a+gc*sqrt(a'*Sigma*a)
DQ_an=Mu+gc*Sigma*a/sqrt(a'*Sigma*a);
ContrQ_an=a.*DQ_an;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure
subplot(2,1,1)
bar(ContrQ_an)
xlim([0 N+1])
title('contributions to VaR, analytical')

subplot(2,1,2)
bar(ContrQ_simul)
xlim([0 N+1])
title('contributions to VaR, simulations')