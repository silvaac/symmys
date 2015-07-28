% this script implements the Expectation-Maximization (EM) algoritm, which estimates 
% the parameters of a multivariate normal distribution when some observations are randomly missing
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input data
load ('db_SwapRatesUS');
Y=[Data(:,2) Data(:,3) Data(:,4)];

DropNum=40; % number of observations dropped

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Original_Series=Y(2:end,:)-Y(1:end-1,:);
[T,N]=size(Original_Series);

% drop data at random
Series_Missing=reshape(Original_Series,T*N,1);
DropIndex=unique(ceil(T*N*rand(DropNum,1)));
Series_Missing(DropIndex)=NaN;
Series_Missing=reshape(Series_Missing,T,N);

% run EM
[E_EM, S_EM , Recovered_Series]=EMestimator(Series_Missing);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display results
Dates=Data(2:end,1);
figure
for n=1:N
    Drop=isnan(Series_Missing(:,n));
    Bad_Dates=Dates(Drop);
    
    Keep=~isnan(Series_Missing(:,n));
    Good_Dates=Dates(Keep);
    
    subplot(N,1,n)
    h1=plot(Good_Dates,Original_Series(Keep,n),'.');
    hold on
    h2=plot(Bad_Dates,Original_Series(Drop,n),'.');
    set(h2,'color','g','markersize',15);
    hold on
    h3=plot(Bad_Dates,Recovered_Series(Drop,n),'.');
    set(h3,'color','r','markersize',14);

    xlim([Dates(1) Dates(end)])
    datetick('x','mmmyy','keeplimits','keepticks');
    grid on
end
legend([h2 h3],'original data','EM-recovered','location',[0.2 .03 .05 .03])