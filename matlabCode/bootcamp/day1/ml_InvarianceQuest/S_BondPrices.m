% this script shows that zero-bond prices cannot be market invariants because 
% their time series converges to a deterministic value at maturity
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; close all; clear;

Start='15-Aug-1999';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis with specific bonds

load('DB_FixedIncome'); % loads structures of specific bonds:  Aug_15_2002, Aug_31_2002 and Nov_15_2002
                        %     with fields:  Expiry: (e.g. 731443), Dates (vector of double) and Yields (vector of double)

Times=find(Aug_15_2002.Dates  >= datenum(Start));
Aug_15_2002.Dates=Aug_15_2002.Dates(Times);
Aug_15_2002.Yields=Aug_15_2002.Yields(Times);
% compute price of zero-coupon equivalent
Aug_15_2002.ZeroPrices = exp(-Aug_15_2002.Yields.*(Aug_15_2002.Expiry-Aug_15_2002.Dates)/365);

Times=find(Nov_15_2002.Dates  >= datenum(Start));
Nov_15_2002.Dates=Nov_15_2002.Dates(Times);
Nov_15_2002.Yields=Nov_15_2002.Yields(Times);
% compute price of zero-coupon equivalent
Nov_15_2002.ZeroPrices = exp(-Nov_15_2002.Yields.*(Nov_15_2002.Expiry-Nov_15_2002.Dates)/365);

% plot series
figure
h=plot(Aug_15_2002.Dates,Aug_15_2002.ZeroPrices);
hold on
h=plot(Aug_15_2002.Expiry,1,'.');
hold on
h=plot(Nov_15_2002.Dates,Nov_15_2002.ZeroPrices);
hold on
h=plot(Nov_15_2002.Expiry,1,'.');
datetick('x','mmmyy','keeplimits','keepticks');
grid on