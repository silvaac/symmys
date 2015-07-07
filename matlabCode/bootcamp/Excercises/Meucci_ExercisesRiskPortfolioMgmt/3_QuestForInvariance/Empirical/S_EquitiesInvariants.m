% this script performs the quest for invariance in the stock market
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: pick one stock from database
Stock_Index=20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('DB_Equities');    % load daily stock prices from the utility sector in the S&P 500
P=Prices(632:end,Stock_Index); % select data after 1/8 rule

% quest for invariance
X=P(2:end)./P(1:end-1);
IIDAnalysis(X);
axisLimits = axis;
text(axisLimits(1:2)*[-0.1,1.1]', axisLimits(3:4)*[0.1,0.9]', ...
    'Analysis for X');

Y=P(2:end)-P(1:end-1);
IIDAnalysis(Y);
axisLimits = axis;
text(axisLimits(1:2)*[-0.1,1.1]', axisLimits(3:4)*[0.1,0.9]', ...
    'Analysis for Y');

Z=X.^2;
IIDAnalysis(Z);
axisLimits = axis;
text(axisLimits(1:2)*[-0.1,1.1]', axisLimits(3:4)*[0.1,0.9]', ...
    'Analysis for Z');

W=P(3:end)-2*P(2:end-1)+P(1:end-2);
IIDAnalysis(W);
axisLimits = axis;
text(axisLimits(1:2)*[-0.1,1.1]', axisLimits(3:4)*[0.1,0.9]', ...
    'Analysis for W');

