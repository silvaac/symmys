% this script shows how to switch between 
% - return-based / relative portfolio weight representation of the efficient frontier
% - price-based / number of securities representation of the efficient frontier
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear; clc; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Horizon=1;

N=30;
Rho=.7;
Current_Prices=rand(N,1);
sigma_Low =.05; sigma_High =.4; 
Step=(sigma_High-sigma_Low)/(N-1); 
sigma=[sigma_Low : Step : sigma_High]';
Mu=.5*sigma;

Budget=1000;

NumPortf=50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Correlation = (1-Rho) * eye(N) + Rho * ones(N,N);
Sigma = diag(sigma)*Correlation*diag(sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[LinRets_ExpectedValues,LinRets_Covariance]=Log2Lin(Horizon*Mu,Horizon*Sigma);
[ExpectedReturn,StandardDeviation, Composition1_Weights]=EfficientFrontier(NumPortf, LinRets_Covariance, LinRets_ExpectedValues);
ExpectedPrice=Budget*(1+ExpectedReturn);
PriceStandardDeviation=Budget*StandardDeviation;
Composition1_Shares=Budget*Composition1_Weights./(ones(NumPortf,1)*Current_Prices');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot prices frontier and absolute composition 
figure

subplot(2,1,1)
h=plot(PriceStandardDeviation,ExpectedPrice);
set(h,'linewidth',2,'color','k');
set(gca,'xlim',[PriceStandardDeviation(1) PriceStandardDeviation(end)]);
grid on
xlabel('standard deviation of portfolio value')
ylabel('expected portfolio value')

subplot(2,1,2)
Data=cumsum(Composition1_Shares,2);
for n=1:N
    x=[PriceStandardDeviation(1); PriceStandardDeviation; PriceStandardDeviation(end)];
    y=[0; Data(:,N-n+1); 0];
    hold on
    h=fill(x,y,[.85 .85 .85]-mod(n,2)*[.2 .2 .2]);
end
set(gca,'xlim',[PriceStandardDeviation(1) PriceStandardDeviation(end)],'ylim',[0 max(max(Data))])
xlabel('standard deviation of portfolio value')
ylabel('number of securities')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot returns frontier and relative composition 
figure

subplot(2,1,1)
h=plot(StandardDeviation,ExpectedReturn);
set(h,'linewidth',2,'color','k');
set(gca,'xlim',[StandardDeviation(1) StandardDeviation(end)]);
grid on
xlabel('standard deviation of portfolio return')
ylabel('expected portfolio return')

subplot(2,1,2)
Data=cumsum(Composition1_Weights,2);
for n=1:N
    x=[StandardDeviation(1); StandardDeviation; StandardDeviation(end)];
    y=[0; Data(:,N-n+1); 0];
    hold on
    h=fill(x,y,[.85 .85 .85]-mod(n-1,2)*[.2 .2 .2]);
end
set(gca,'xlim',[StandardDeviation(1) StandardDeviation(end)],'ylim',[0 max(max(Data))])
xlabel('standard deviation of portfolio return')
ylabel('relative weight')