% since the certainty equivalent for a power utility functions is homogeneous
% the Euler decomposition of satisfaction into the contributions from each
% exposure applies. The contributions are computed both as expecations and as empirical gradient
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs (log-normal market, power utility)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of assets
N=40; 

% log-normal copula 
r=.99; 

% number of simulations
J=10000;

% risk aversion
g=.3;

% allocation
a=rand(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulate market at the investment horizon
m=.1*rand(N,1);
s=.2*rand(N,1);
Rho=(1-r)*eye(N)+r*ones(N,N);
S=diag(s)*Rho*diag(s);
% compounded returns
C=mvnrnd(m,S,J);
% prices
P=exp(C);

% generate sample of final wealth and tweak values for empirical derivative
W=P*a;
% compute expected utility and certainty equivalent
u=mean(W.^g);
CE=u^(1/g);
% compute tweaked expected utility and certainty equivalent for gradient
D_CE=[];
for n=1:N
    e=zeros(N,1); 
    e(n)=.01;
    W_up=P*(a+e);    
    u_up=mean(W_up.^g);
    CE_up=u_up^(1/g);

    D_CE=[D_CE
        (CE_up-CE)/e(n)];
end
% compute marginal contributions
Contr=a.*D_CE

% compute gradient with expectation formula
D_CE_exp=[];
for n=1:N
    D_CEn=(CE^(g-1))*mean(P(:,n).*(W.^(1-g)));
    D_CE_exp=[D_CE_exp
        D_CEn];
end
% compute marginal contributions with expectation formula
Contr_exp=a.*D_CE_exp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,1,1)
bar(Contr_exp)
xlim([0 N+1])
grid on
title('marginal contributions: expectation formula')

subplot(2,1,2)
bar(Contr)
xlim([0 N+1])
grid on
title('marginal contributions: empirical gradient')
