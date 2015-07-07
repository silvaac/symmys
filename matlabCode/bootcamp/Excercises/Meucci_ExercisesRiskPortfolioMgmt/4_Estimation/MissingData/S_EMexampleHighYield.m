% this script implements the Expectation-Maximization (EM) algoritm, which estimates 
% the parameters of a multivariate normal distribution when some observations are randomly missing
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input data
load ('db_HighYieldIndices');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T,N]=size(Data);
Series=log(Data(2:end,:))-log(Data(1:end-1,:));
NANs_Index=find(abs(Series)>10^10);
Series(NANs_Index)=NaN;

% run EM
[E_EM, S_EM , Recovered_Series]=EM(Series);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display results

figure
for n=1:N
    Drop=isnan(Series(:,n));
    Bad_Dates=Dates(Drop);
    
    Keep=~isnan(Series(:,n));
    Good_Dates=Dates(Keep);
    
    subplot(N,1,n)
    h1=plot(Good_Dates,Series(Keep,n),'.');
    hold on
    h3=plot(Bad_Dates,Recovered_Series(Drop,n),'.');
    set(h3,'color','r','markersize',14);

    xlim([Dates(1) Dates(end)])
    datetick('x','mmmyy','keeplimits','keepticks');
    grid on
end
legend(h3,'EM-recovered data','location',[0.2 .03 .05 .03])