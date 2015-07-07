% this script performs the quest for invariance in the derivatives market
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; close all; clear;
% load implied vol for options on SPX for different time to maturity and moneyness
load('DB_Derivatives');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple univariate test

% select the implied vol time series for a specific time to maturity and moneyness 
maturityIndex=1;    % 1..6
moneynessIndex=4;   % 1..7

% quest for invariance
X=diff(impVol(1:5:end,maturityIndex,moneynessIndex));
IIDAnalysis(X);
axisLimits = axis;
text(axisLimits(1:2)*[-0.1,1.1]', axisLimits(3:4)*[0.1,0.9]', 'Changes in implied vol');

Y=diff(log(impVol(1:5:end,maturityIndex,moneynessIndex)));
IIDAnalysis(Y);
axisLimits = axis;
text(axisLimits(1:2)*[-0.1,1.1]', axisLimits(3:4)*[0.1,0.9]', 'Changes in log of implied vol');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multivariate test with AR(1) structure

[T,Mat,Mon]=size(impVol(1:5:end,:,:));
Z=reshape(log(impVol(1:5:end,:,:)),T,Mat*Mon);
X=Z(2:end,:);
F=[ones(T-1,1) Z(1:end-1,:)];
E_XF=X'*F/T;
E_FF=F'*F/T;
B=E_XF*inv(E_FF);
Eps=X-F*B';

IIDAnalysis(Eps(:,3));
axisLimits = axis;
text(axisLimits(1:2)*[-0.1,1.1]', axisLimits(3:4)*[0.1,0.9]', 'VAR(1) residuals');
