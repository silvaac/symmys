% this script plots the joint and marginal distributions of the swap rate daily changes
% it computes the summary statistics of the distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs

% upload sample
load('db_SwapRatesUS');

% select invariants (two among USD SWAP 5Y rate, USD SWAP 5Y rate, USD SWAP 5Y rate, daily changes)
Pick=[1 2]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dates=Data(:,1);
% plots series
figure
plot(Dates,Data(:,2),'r');
hold on
plot(Dates,Data(:,3),'g');
hold on
plot(Dates,Data(:,4),'b');
legend(Fields(2).Name,Fields(3).Name,Fields(4).Name,'location','best')
xlim([Dates(1) Dates(end)])
datetick('x','mmmyy','keeplimits','keepticks');
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = Data(2:end,[2:end])-Data(1:end-1,[2:end]);
T=size(X,1);


figure 
% marginals
NumBins=round(10*log(T));

subplot('Position',[.05 .3 .2 .6]) 
[n,D]=hist(X(:,2),NumBins);
barh(D,n,1);
[y_lim]=get(gca,'ylim');
set(gca,'xtick',[])
grid on

subplot('Position',[.3 .05 .6 .2]) 
[n,D]=hist(X(:,1),NumBins);
bar(D,n,1);
[x_lim]=get(gca,'xlim');
set(gca,'ytick',[])
grid on

% scatter plot
subplot('Position',[.3 .3 .6 .6]) 
h=plot(X(:,1),X(:,2),'.');
set(gca,'xlim',x_lim,'ylim',y_lim)
grid on
xlabel(Fields(Pick(1)+1).Name)
ylabel(Fields(Pick(2)+1).Name)

% compute summary statistics

Sample_Mean=mean(X)
Sample_Covariance=cov(X);
Sample_Standard_Deviation=sqrt(diag(Sample_Covariance))'
Sample_Correlation=Sample_Covariance(1,2)/(sqrt(Sample_Covariance(1,1)*Sample_Covariance(2,2)))
