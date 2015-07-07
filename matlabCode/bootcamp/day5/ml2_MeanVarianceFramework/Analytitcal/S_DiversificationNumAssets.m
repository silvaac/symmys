% this script shows the diversification effect of increasing the number of securities in a market
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_Low=0.05;
E_High=0.10;
S_Low=0.10;
S_High=0.25;
NumAssets=[2 5 10 20 100];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
for k=1:length(NumAssets)
    StepE=(E_High-E_Low)/(NumAssets(k)-1);
    StepS=(S_High-S_Low)/(NumAssets(k)-1);
    E=[E_Low:StepE:E_High]';
    s=[S_Low:StepS:S_High];
    S=diag(s.^2);
    NumPortf=20;
    [ExpectedReturn,Volatility, Composition]=EfficientFrontier(NumPortf, S, E);
    plot(Volatility, ExpectedReturn,'color','k');
    axis([0 S_High*1.1 E_Low E_High*1.1]);
    hold on
end
grid on

