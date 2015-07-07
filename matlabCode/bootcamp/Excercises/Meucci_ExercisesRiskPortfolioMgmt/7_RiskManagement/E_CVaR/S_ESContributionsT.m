% this script computes the expected shortfall and the contributions to ES from each security 
% - analytically, under the Student t assumption for the market
% - in simulations, using the conditional expectation definition of the contributions
% See Sec 5.6 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs 

% number of assets
N=40; 

% market parameters (student t distribution)
Mu=rand(N,1);
A=rand(N,N)-.5;
Sigma=A*A';
nu=7;

% allocation
a=rand(N,1)-.5;

% ES confidence
c=.95;

NumSimul=100000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate market scenarios
l=ones(NumSimul,1);
diag_s=diag(sqrt(diag(Sigma)));
C=inv(diag_s)*Sigma*inv(diag_s);
X=mvtrnd(C,nu,NumSimul/2);
X=[X         % symmetrize simulations
    -X];
M = l*Mu' + X*diag_s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulations

% compute the objective
Psi = M*a;

% compute cut-off spectrum (step function) for empirical ES estimation, see (5.218)
th=ceil((1-c)*NumSimul); % threshold
spc=zeros(NumSimul,1);
spc([1:th])=1;
spc=spc/sum(spc);

% compute ES
[Sort_Psi,Index] = sort(Psi);
ES_simul=Sort_Psi'*spc

% sort market according to order induced by objective's realizations
Sort_M=M(Index,:);
for n=1:N
    % compute gradient as conditional expectation
    DES_simul(n)=spc'*Sort_M(:,n);
end
% compute contributions
ContrES_simul=a.*DES_simul';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical

% this does NOT depend on the allocation...
ES_standardized=1/(1-c)*quad('tinv',10^(-8),1-c,[],[],nu); 

% ...the dependence on the allocation is analytical
ES_an=Mu'*a+ES_standardized*sqrt(a'*Sigma*a)
DES_an=Mu+ES_standardized*Sigma*a/sqrt(a'*Sigma*a);
ContrES_an=a.*DES_an;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

figure

subplot(2,1,1)
bar(ContrES_an)
xlim([0 N+1])
title('contributions to ES, analytical')

subplot(2,1,2)
bar(ContrES_simul)
xlim([0 N+1])
title('contributions to ES, simulations')