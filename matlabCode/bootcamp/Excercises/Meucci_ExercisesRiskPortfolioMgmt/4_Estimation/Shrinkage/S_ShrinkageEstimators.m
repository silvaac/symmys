% this script computes the multivariate shrinkage estimators of location and scatter under the normal assumption
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
N=5;
T=100;

Mu=rand(N,1);
A=rand(N,N)-.5;
Sigma=A*A';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate normal sample
X=mvnrnd(Mu,Sigma,T);

% estimate sample parameters
Mu_hat=mean(X)';
Sigma_hat=cov(X)*(T-1)/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shrinkage of location parameter 

% target 
b=0*ones(N,1); 

% compute optimal weight
Lambda_hat=eig(Sigma_hat); 
a=1/T*(sum(Lambda_hat)-2*max(Lambda_hat))/((Mu_hat-b)'*(Mu_hat-b)); 

% restrict to sensible weight
a=max(0,min(a,1));     

% shrink
Mu_shr=(1-a)*Mu_hat+a*b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shrinkage of scatter parameter

% target
C=mean(Lambda_hat)*eye(N);

% compute optimal weight
Numerator=0;
for t=1:T
    Numerator = Numerator +  1/T*trace(  (X(t,:)'*X(t,:)-Sigma_hat)^2  ) ;
end
Denominator=trace( (Sigma_hat-C)^2);
a=1/T*Numerator/Denominator;

% restrict to sensible weight
a=max(0,min(a,1)); 

% shrink
Sigma_shr=(1-a)*Sigma_hat+a*C;