% this script display the empirical copula of a set of market variables
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters 
load db_FX 
Display=[1 2];  % 1 = Spot USD/EUR; 2 = Spot USD/GBP; 3 = Spot USD/JPY; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=log(Data(2:end,2:end))-log(Data(1:end-1,2:end));

% compute empirical copula by sorting
[NumObs,K]=size(X);
[W,C]=sort(X);
for k=1:K
    x=C(:,k);
    y=[1:NumObs];
    xi=[1:NumObs];
    yi = interp1(x,y,xi);
    Copula(:,k)=yi/(NumObs+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
% marginals
NumBins=round(10*log(NumObs));

subplot('Position',[.05 .3 .2 .6]) 
[n,D]=hist(X(:,Display(2)),NumBins);
barh(D,n,1);
grid on

subplot('Position',[.3 .05 .6 .2]) 
[n,D]=hist(X(:,Display(1)),NumBins);
bar(D,n,1);
grid on

% scatter plot
subplot('Position',[.3 .3 .6 .6]) 
h=plot(Copula(:,Display(1)),Copula(:,Display(2)),'.');
grid on
title('Copula')
xlabel(num2str(Fields(Display(1)+1).Name));
ylabel(num2str(Fields(Display(2)+1).Name));