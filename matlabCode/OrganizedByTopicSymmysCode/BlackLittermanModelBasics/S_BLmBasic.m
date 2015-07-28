% this script describes to basic market-based Black-Litterman approach 
% in particular: full confidence = conditional
%                no confidence = reference model
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clear; clc; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pick = [1; 6];              % views pick matrix
NumPortf=40;                   % number of MV-efficient portfolios 
load('CovNRets');               % input Covariance and Mu of asset returns from database...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumAssets=size(Sigma,2);
NumViews=length(Pick);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute efficient frontier
[M,S]=Log2Lin(Mu,Sigma);
[ExpectedReturn,Volatility, Composition]=EfficientFrontier(NumPortf, S, M);
PlotFrontier(Composition)

% modify expected returns the Black-Litterman way and compute new efficient frontier 
P = [1 0 0 0 0 -1];
Omega=P*Sigma*P';
Views = sqrt(diag(Omega));   % views value

[BLMu,BLSigma]=BLmFormulas(Mu,Sigma,P,Views,Omega)

[M,S]=Log2Lin(BLMu,BLSigma);                
[E,V, BLComposition]=EfficientFrontier(NumPortf, S, M);
PlotFrontier(BLComposition)