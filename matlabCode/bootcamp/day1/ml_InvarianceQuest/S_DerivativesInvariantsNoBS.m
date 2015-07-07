% Invariants for options
clear; clc; close all

% load historical price series for a call option on SPX
% maturity: June 2009.
load SPQ_GT_prices dates prices spx_prices;

% IID analysis for call prices
IIDAnalysis(dates,prices)

% IID analysis for changes in call price
priceChange = prices(2:end) - prices(1:end-1);
%IIDAnalysis2(priceChange)

% IID analysis for total returns
Ret = prices(2:end)./prices(1:end-1);
%IIDAnalysis2(Ret)

% load implied vol historicals for different moneyness Delta=K/spot and
% different times to maturity
clear
r_free = 0.04;
load HistoricalVolSurface;

% convert implied vol to prices of call options
callPrice = zeros(length(spot),length(Delta));
for i = 1:length(Delta),
    callPrice(:,i) = BlackScholesCall(spot,Delta(i)*spot,r_free,impVol(:,i),days2Maturity(i)/252);
end

% Compute historical pseudo-returns for different Delta and time to maturity
Z = callPrice(2:end,:)./callPrice(1:end-1,:);

% IID analysis for pseudo-returns
optionIndex = 17; % 1..42
IIDAnalysis2(Z(:,optionIndex));
