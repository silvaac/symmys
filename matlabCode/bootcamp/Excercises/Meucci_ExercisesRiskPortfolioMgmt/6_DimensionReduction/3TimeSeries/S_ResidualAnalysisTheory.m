% Script for explicit factors and implicit loadings: 
% Analysis of residual
clear; clc; close all

%------------------------------
% Input parameters:
N = 3;    % number of stocks
K = 2;     % number of factors
tau = 12/12; % years
numSimulations = 100000;
useNonRandomFactor = 1
useZeroMeanLinearFactors = 0
%------------------------------

% generate covariance matrix with volatility 25%
Diag=.25*ones(N+K,1);
if useNonRandomFactor,
    Diag(end)=10^(-26); % make one factor deterministic -> add constant
end
dd=rand(N+K,N+K)-.5;
dd=dd*dd';
[dd,C]=cov2corr(dd);

covJointXF = diag(Diag)*C*diag(Diag);

SigmaX = covJointXF(1:N,1:N);
SigmaXF = covJointXF(1:N,(N+1):(N+K));
SigmaF = covJointXF((N+1):(N+K),(N+1):(N+K));

% generate mean returns for stock and factors ~ 10% per annum
muX = 0.1*randn(N,1); 
muF = 0.1*randn(K,1); 

% statitics of linear returns analytically, since Y = 1+[R; Z] is lognormal
mu = [muX; muF];
E_Y = exp(mu*tau+diag(covJointXF*tau)/2);
E_YY = (E_Y*E_Y').*exp(covJointXF*tau);
E_R = E_Y(1:N) - ones(N,1);
E_Z = E_Y((N+1):end) - ones(K,1);
E_RR = E_YY(1:N,1:N) - ones(N,1)*E_R' - E_R*ones(1,N) - ones(N,N);
E_ZZ = E_YY((N+1):end,(N+1):end) - ones(K,1)*E_Z' - E_Z*ones(1,K) - ones(K,K);
E_RZ = E_YY(1:N,(N+1):end) - ones(N,1)*E_Z' - E_R*ones(1,K) - ones(N,K);
SigmaZ = E_ZZ - E_Z*E_Z';
SigmaR = E_RR - E_R*E_R';
SigmaRZ = E_RZ - E_R*E_Z';

% generate Monte Carlo simulations
sims = mvnrnd(mu*tau, covJointXF*tau, numSimulations);
X = sims(:,1:N);
F = sims(:,(N+1):end);
R = exp(X)-1;
Z = exp(F)-1;
if useZeroMeanLinearFactors,
    % enforce Z sample to be zero-mean, equivalent to muF = -diag(SigmaF)/2
    Z = Z - repmat(mean(Z),numSimulations,1);
end

% compute sample estimates
E_R_hat = mean(R,1)';
E_Z_hat = mean(Z,1)';
SigmaR_hat = cov(R,1);
SigmaZ_hat = cov(Z,1);

% compute simulation errors
errSigmaR = norm(SigmaR-SigmaR_hat,'fro')/norm(SigmaR,'fro');
fprintf('Simulation error in sample cov(R) as a percentage on true cov(R) = %.1f%%\n',errSigmaR*100);
errSigmaZ = norm(SigmaZ-SigmaZ_hat,'fro')/norm(SigmaZ,'fro');
fprintf('Simulation error in sample cov(Z) as a percentage on true cov(Z) = %.1f%%\n',errSigmaZ*100);

% compute OLS loadings for the linear return model
B = E_RZ*inv(E_ZZ);
B_hat = R'*Z*inv(Z'*Z);
errB = norm(B-B_hat,'fro')/norm(B,'fro');
fprintf('Simulation error in sample OLS loadings as a percentage on true OLS loadings = %.1f%%\n',errB*100);

U = R - Z*B_hat';
Corr=corr([U Z]);
%Corr=cov([U Z],1);

Corr_U=Corr(1:N,1:N)
Corr_UZ=Corr(1:N,N+1:N+K)

SigmaU_hat = cov(U,1);
BSBplusSu=B_hat*SigmaZ_hat*B_hat'+SigmaU_hat;
errSigmaR_model1 = norm(SigmaR_hat-BSBplusSu,'fro') ...
    /norm(SigmaR_hat,'fro');
fprintf('Sigma_R - (B*Sigma_Z*B + Sigma_U) = %.1f%%\n',errSigmaR_model1*100);