% this script performs the quest for invariance in the fixed income market
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: pick time-to-maturity for one point on the yield curve
ycMaturityIndex=4;       % 1..6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load government yield curve and bond yield data for different dates
load('DB_FixedIncome');

% select the yield time series for a constant time-to-maturity
yield=ycYieldPercent(:,ycMaturityIndex); 

% quest for invariance:
% changes in the yield curve
X = yield(2:end) - yield(1:end-1);
IIDAnalysis(X);
axisLimits = axis;
text(axisLimits(1:2)*[-0.1,1.1]', axisLimits(3:4)*[0.1,0.9]', ...
    'Changes in yield curve');

% changes in the logarithm of the yield curve
Y = log(yield(2:end)) - log(yield(1:end-1));
IIDAnalysis(Y);
axisLimits = axis;
text(axisLimits(1:2)*[-0.1,1.1]', axisLimits(3:4)*[0.1,0.9]', ...
    'Changes in log of yield curve');
