% this script shows how the efficient frontier depends on the investment horizon
% and highlights the problems arising from mistakenly using compounded
% returns instead of linear returns for this analysis
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clear; clc; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Horizons=[1/12 3];

N=8;
Rho=.7;
Current_Prices=rand(N,1);
sigma_Low =.05; sigma_High =.4; 
Step=(sigma_High-sigma_Low)/(N-1);
sigma=[sigma_Low : Step : sigma_High]';
Mu=.5*sigma;

Budget=1000;

NumPortf=20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Correlation = (1-Rho) * eye(N) + Rho * ones(N,N);
Sigma = diag(sigma)*Correlation*diag(sigma);


% shorter horizon analysis
[LinRets_ExpectedValues,LinRets_Covariance]=Log2Lin(Horizons(1)*Mu,Horizons(1)*Sigma);
[ExpectedReturn,StandardDeviation, Composition1_Weights]=EfficientFrontier(NumPortf, LinRets_Covariance, LinRets_ExpectedValues);

% longer horizon analysis
[LinRets_ExpectedValues,LinRets_Covariance]=Log2Lin(Horizons(2)*Mu,Horizons(2)*Sigma);
[ExpectedReturn,StandardDeviation, Composition2_Weights]=EfficientFrontier(NumPortf, LinRets_Covariance, LinRets_ExpectedValues);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,1,1)
Data=cumsum(Composition1_Weights,2);
for n=1:N
    x=[StandardDeviation(1); StandardDeviation; StandardDeviation(end)];
    y=[0; Data(:,N-n+1); 0];
    hold on
    h=fill(x,y,[.85 .85 .85]-mod(n-1,2)*[.2 .2 .2]);
end
axis([StandardDeviation(1) StandardDeviation(end) 0 max(max(Data))]);
xlabel('standard deviation of portfolio return')
ylabel('relative weight')
title(['Investment horizon = ' num2str(Horizons(1)) ' years'],'fontweight','bold');

subplot(2,1,2)
Data=cumsum(Composition2_Weights,2);
for n=1:N
    x=[StandardDeviation(1); StandardDeviation; StandardDeviation(end)];
    y=[0; Data(:,N-n+1); 0];
    hold on
    h=fill(x,y,[.85 .85 .85]-mod(n-1,2)*[.2 .2 .2]);
end
axis([StandardDeviation(1) StandardDeviation(end) 0 max(max(Data))]);
xlabel('standard deviation of portfolio return')
ylabel('relative weight')
title(['Investment horizon = ' num2str(Horizons(2)) ' years'],'fontweight','bold');