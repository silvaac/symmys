% this script describes to basic market-based Black-Litterman approach 
% in particular: full confidence = conditional
%                no confidence = reference model
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clear; clc; close all

load('CovNRets'); % input Covariance and Mu of asset returns from database...

NumAssets=size(Sigma,2);

%% compute efficient frontier

[M,S]=Log2Lin(Mu,Sigma);
NumPortf=40;% number of MV-efficient portfolios 
Inputs.Exp=M;
Inputs.Cov=S;

Constr.Aeq=ones(1,NumAssets);
Constr.beq=1;
Constr.Aleq=-eye(NumAssets);
Constr.bleq=zeros(NumAssets,1);

[ExpectedReturn,Volatility, Composition]=EfficientFrontier(NumPortf, Inputs, Constr);
PlotFrontier(Composition, Volatility)

%% modify expected returns the Black-Litterman way and compute new efficient frontier 
% Pick = [1; 6];              % views pick matrix
% NumViews=length(Pick);
 
P = [1 0 0 0 0 -1];
Omega=P*Sigma*P'; %c=1/2
Views = sqrt(diag(Omega));   % views value

[BLMu,BLSigma]=BLmFormulas(Mu,Sigma,P,Views,Omega);

[InputsBL.Exp,InputsBL.Cov]=Log2Lin(BLMu,BLSigma);                
[E,V, BLComposition]=EfficientFrontier(NumPortf, InputsBL, Constr);
PlotFrontier(BLComposition,V)