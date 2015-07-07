% this script plots the joint and marginal distributions of the swap rate daily changes
% it computes the summary statistics of the distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% upload sample
load('DBUsSwapRates');

% plots series
figure
plot(Time,Series(1).Data,'r');
hold on
plot(Time,Series(2).Data,'k');
hold on
plot(Time,Series(3).Data);
legend(cell2mat(Series(1).Name),cell2mat(Series(2).Name),cell2mat(Series(3).Name))
xlim([Time(1) Time(end)])
datetick('x','mmmyy','keeplimits','keepticks');
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select invariants (two among 2yr, 5yr, 10yr swap rates daily changes)
Pick=[1 3]; 

Rates=[Series(Pick(1)).Data Series(Pick(2)).Data];
X = Rates(2:end,:)-Rates(1:end-1,:);
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
xlabel(Series(Pick(1)).Name)
ylabel(Series(Pick(2)).Name)

% compute summary statistics

Sample_Mean=mean(X)
Sample_Covariance=cov(X);
Sample_Standard_Deviation=sqrt(diag(Sample_Covariance))'
Sample_Correlation=Sample_Covariance(1,2)/(sqrt(Sample_Covariance(1,1)*Sample_Covariance(2,2)))
