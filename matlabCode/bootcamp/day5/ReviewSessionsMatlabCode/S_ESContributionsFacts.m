% this script computes the expected shortfall and the contributions to ES from each factor 
% in simulations, using the conditional expectation definition of the contributions
% See Sec 5.6 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;
% inputs 

K=10; % number of factors
N=30; % number of securities
a=rand(N,1); % allocation
c=.95; % ES confidence

% generate market simulations

% factor loadings 
B=rand(N,K); 

% student t parameters
nu=10;
eps=.1;
g=.1;

A=rand(K,K)-.5;
C=A*A';
[dd ,sigma_f]=cov2corr(C);
sigma_u=toeplitz(exp(-g*(0:29)));
sigma=blkdiag(eps*sigma_f,eps^2*sigma_u);
[diag_sigma, corr]=cov2corr(sigma);

%% scenarios
NumSimul=10000;
X=mvtrnd(corr,nu,NumSimul);
X=X*diag(diag_sigma);
X=exp(X);
F=X(:,1:K);
U=X(:,K+1:end); 
U=U-repmat(mean(U),NumSimul,1);
M=F*B'+U;

%% risk management

% compute the objective
Psi = M*a; 

% compute ES
th=ceil((1-c)*NumSimul); % threshold
spc=zeros(NumSimul,1);
spc([1:th])=1;
spc=spc/sum(spc);

[Sort_Psi,Index] = sort(Psi);
ES_simul=Sort_Psi'*spc;

%or: ES_simul=mean(Psi(Psi<=quantile(Psi,1-c)));

% augment factor set to include residual
F_=[F U*a];
% compute portfolio-level loadings
b_=[a'*B 1];

% compute contributions
ContrES_simul=b_.*mean(F_(Psi<=quantile(Psi,1-c),:));

%% plots
figure
bar(ContrES_simul)
xlim([0 K+2])
title('contributions to ES')