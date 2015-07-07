% this script computes the expected shortfall and the contributions to ES from each factor 
% in simulations, using the conditional expectation definition of the contributions
% See Sec 5.6 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs 

K=10; % number of factors
N=40; % number of securities
a=rand(N,1); % allocation
c=.95; % ES confidence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate market simulations

% factor loadings 
B=rand(N,K); 

% student t parameters
nu=15;
eps=.1;
g=.1;

A=rand(K,K)-.5;
C=A*A';
[dd ,sigma_f]=cov2corr(C);
sigma_u=eye(N);
for n=1:N-1
    sigma_u=sigma_u+exp(-g*n)*(diag(ones(N-n,1),n)+diag(ones(N-n,1),-n));
end
sigma=blkdiag(eps*sigma_f,eps^2*sigma_u);
[diag_sigma, corr]=cov2corr(sigma);

% scenarios
NumSimul=10000;
l=ones(NumSimul,1);
X=mvtrnd(corr,nu,NumSimul/2);
X=[X         % symmetrize simulations
    -X];
X=X*diag(diag_sigma);
X=exp(X);
F=X(:,1:K);
U=X(:,K+1:end); 
U=U-l*mean(U);
M=F*B'+U;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% risk management

% compute the objective
Psi = M*a; 

% compute ES
th=ceil((1-c)*NumSimul); % threshold
spc=zeros(NumSimul,1);
spc([1:th])=1;
spc=spc/sum(spc);

[Sort_Psi,Index] = sort(Psi);
ES_simul=Sort_Psi'*spc;

% augment factor set to include residual
F_=[F U*a];
% compute portfolio-level loadings
b_=[a'*B 1];
% sort factors according to order induced by objective's realizations
Sort_F_=F_(Index,:);
for k=1:K+1
    % compute gradient as conditional expectation
    DES_simul(k)=spc'*Sort_F_(:,k);
end
% compute contributions
ContrES_simul=b_.*DES_simul;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure
bar(ContrES_simul)
xlim([0 K+2])
title('contributions to ES')