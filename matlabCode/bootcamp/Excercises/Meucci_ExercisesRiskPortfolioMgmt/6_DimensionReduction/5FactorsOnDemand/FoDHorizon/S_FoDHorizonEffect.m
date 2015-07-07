% Script to study horizon effect on Factors on Demand approach to beta computation
% see Meucci, A. "Common Misconceptions About 'Beta' - Hedging, Estimation and Horizon Effects" 
% GARP's Risk Professional Magazine, June 2010 
% available at http://ssrn.com/abstract=1619923


clear; clc; close all

% Compounded returns follow the linear model X = tau*muX + D*F + epsilon, where 
% tau: investment horizon (in weeks)
% muX: expected weekly compounded returns
% F: factor compounded returns, with zero expectation and tau-proportional covariance
% D: matrix of factor loadings
% epsilon: uncorrelated (idiosyncratic) shocks.
% R=exp(X)-1 and Z=exp(F)-1 are the linear returns

% Load parameters of the model: D, muX, sigmaF, sigmaEps
load db_LinearModel;

% Specify range of investment horizon, weeks
tauRangeWeeks = 1:52;

% constants
[N,K] = size(D);

aMinusTauMuX = zeros(1,length(tauRangeWeeks));
normDminusB = zeros(1,length(tauRangeWeeks));
minCorrU = zeros(1,length(tauRangeWeeks));
meanCorrU = zeros(1,length(tauRangeWeeks));
maxCorrU = zeros(1,length(tauRangeWeeks));
for i = 1:length(tauRangeWeeks),
    tau = tauRangeWeeks(i);
    
    % statitics of linear returns (analytically)
    % let Y = [1+R; 1+Z] ~ LogN(mu*tau, covJointXF*tau)
    mu = [muX; zeros(K,1)];
    % covariance of [X; F] for tau=1:
    covJointXF = [D*sigmaF*D'+sigmaEps, D*sigmaF
                  sigmaF*D', sigmaF];
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

    % compute OLS loadings for the linear return model
    B = SigmaRZ/SigmaZ;
    a = E_R - B*E_Z;
    aMinusTauMuX(i) = norm(a - tau*muX);
    normDminusB(i) = norm(D-B,'fro');
    
    % pairwise correlations of U
    SigmaU = SigmaR - B*SigmaRZ';
    [stdU,corrU] = cov2corr(SigmaU);
    stackedCorrU = corrU(:);
    minCorrU(i) = min(min(abs(corrU)));
    meanCorrU(i) = (N*mean(mean(abs(corrU)))-1)/(N-1);
    corrU(1:(N+1):N*N) = 0;
    maxCorrU(i) = max(max(abs(corrU)));
end

figure
plot(tauRangeWeeks, aMinusTauMuX);
grid on
xlabel('Investment horizon, \tau, weeks')
title('Norm of (a-\tau\mu_X)')

figure
plot(tauRangeWeeks, normDminusB)
grid on
xlabel('Investment horizon, \tau, weeks')
title('Norm of (D-B)')

figure
plot(tauRangeWeeks, maxCorrU, 'r')
hold on
plot(tauRangeWeeks, meanCorrU, 'b')
plot(tauRangeWeeks, minCorrU, 'g')
title('Pairwise correlations of elements of U')
xlabel('Investment horizon, \tau, weeks')
legend('Max absolute corr','Mean absolute corr','Min absolute corr',0)

