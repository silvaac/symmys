% this script performs the quest for invariance in the bond market
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; close all; clear;

Start='15-Aug-1999';
Stop='04-Jan-2003';
TimeToMaturity=10;   % select rolling-horizon analysis among 2, 5 or 10 years

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis with rolling zero-coupon prices 
[Dates4ZeroPrices,ZeroPrices]=ProcessDB(TimeToMaturity,Start,Stop);
IIDAnalysis(Dates4ZeroPrices,ZeroPrices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis with pseudo-returns of rolling zero-coupon prices 
ZeroReturns=[ZeroPrices(2:end)./ZeroPrices(1:end-1)];
Dates4ZeroReturns=Dates4ZeroPrices(2:end);
IIDAnalysis(Dates4ZeroReturns,ZeroReturns);