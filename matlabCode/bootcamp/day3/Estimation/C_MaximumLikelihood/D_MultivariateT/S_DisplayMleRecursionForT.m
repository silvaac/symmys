% this script demonstrates the recursive ML estimation of the location and scatter
% parameters of a multivariate Student t distribution
% see Section 4.3.1 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci


clear; clc; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
N=4;
T=10000;

nu=1;
Mu=rand(N,1);
A=rand(N,N)-.5;
Sigma=A*A';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=sqrt(diag(Sigma));
C=diag(1./s)*Sigma*diag(1./s);
X=ones(T,1)*Mu'+ mvtrnd(C,nu,T)*diag(s);

Tolerance=10^(-10);
[Mu_hat,Sigma_hat] = MleRecursionForT(X,nu,Tolerance);

Mu
Mu_hat

Sigma
Sigma_hat

