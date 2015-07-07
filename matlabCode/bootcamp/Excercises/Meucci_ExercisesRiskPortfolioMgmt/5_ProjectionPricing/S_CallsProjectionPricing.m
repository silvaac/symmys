% this script projects the distribution of the market invariants for the derivatives market 
% Then it computes the distribution of prices at the investment horizon
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%
% inputs
tau_tilde=5;    % estimation step (days)
tau=40;         % time to horizon (days)
Time2Mats=[100 150 200 250 300];    % current time to maturity of call options in days
Strikes = [850 880 910 940 970];    % strikes of call options, same dimension as Time2Mat

r_free=0.04;    % risk-free rate
J=10000;       % number of simulations

%%%%%%%%%%%%%%%%%%%
% load 'spot' for underlying and current vol surface, given by
% 'impVol' for different 'days2Maturity' and 'moneyness' (K/S)
load DB_ImplVol;
numCalls = length(Time2Mats);
timeLength = length(spot);
numSurfPoints = length(days2Maturity)*length(moneyness);

%%%%%%%%%%%%%%%%%%%
% estimate invariant distribution assuming normality
% variables in X are changes in log(spot) and changes in log(imp.vol)
% evaluated at the 'numSurfPoints' points on the vol surface (vectorized).
X = zeros(timeLength-1,numSurfPoints+1);
% log-changes of underlying spot
X(:,1) = diff(log(spot));

% log-changes of implied vol for different maturities
impVolSeries = reshape(impVol,timeLength,numSurfPoints);
for i = 1:numSurfPoints,
    X(:,i+1) = diff(log(impVolSeries(:,i)));
end
muX = mean(X);
SigmaX = cov(X,1);

%%%%%%%%%%%%%%%%%%%
% project distribution to investment horizon
muX = muX*tau/tau_tilde;
SigmaX = SigmaX*tau/tau_tilde;

%%%%%%%%%%%%%%%%%%%
% linearly interpolate the vol surface at the current time to obtain
% implied vol for the given calls today, and price the calls
spot_T = spot(end);
volSurf_T = squeeze(impVol(end,:,:));
time2Mat_T = Time2Mats;
moneyness_T = Strikes/spot_T;
impVol_T = interpne(volSurf_T,[time2Mat_T',moneyness_T'],{days2Maturity,moneyness})';
callPrice_T = BlackScholesCall(spot_T,Strikes,r_free,impVol_T,Time2Mats/252);

%%%%%%%%%%%%%%%%%%%
% generate simulations at horizon
X_ = mvnrnd(muX,SigmaX,J);

% interpolate vol surface at horizon for the given calls
spot_ = spot_T*exp(X_(:,1));
impVol_ = zeros(J,numCalls);
for j = 1:J
    volSurf = volSurf_T.*exp(reshape(X_(j,2:end),length(days2Maturity),length(moneyness)));
    time2Mat_ = Time2Mats-tau;
    moneyness_ = Strikes/spot_(j);
    impVol_(j,:) = interpne(volSurf,[time2Mat_',moneyness_'],{days2Maturity,moneyness})';
end

% price the calls at the horizon
callPrice_ = zeros(J,numCalls);
for i = 1:numCalls
    callPrice_(:,i) = BlackScholesCall(spot_,Strikes(i),r_free,impVol_(:,i),time2Mat_(i)/252);
end
LinearRets = callPrice_./repmat(callPrice_T,J,1) - 1;

NumBins=round(10*log(J));
for i = 1:numCalls
    figure
    subplot(2,1,1)
    hist(callPrice_(:,i),NumBins)
    xlabel('call price')
    subplot(2,1,2)
    plot(spot_,callPrice_(:,i),'.')
    xlabel('spot price')
    ylabel('call price')
end