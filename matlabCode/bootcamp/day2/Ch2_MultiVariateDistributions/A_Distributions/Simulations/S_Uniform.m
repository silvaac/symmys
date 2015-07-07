    % this script simulates a bivariate uniform distribution on an ellipsoid
% it shows its pdf as proxied by the 3-D histogram
% it shows the distribution of a generic linear combination (portfolio) of the two variables
% it computes the summary statistics of the market
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters

Mu=[0 0]';
s=[0.3 0.2]';
r=.5;

NumSimulations=100000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=[s(1)^2     r*s(1)*s(2)
    r*s(1)*s(2)    s(2)^2];

NumSimulParallelotope=4/pi*prod(sqrt(diag(Sigma)))/sqrt(det(Sigma))*NumSimulations;

% generate uniform sample on the parallelotope
U_1 = Mu(1)+2*s(1)*(rand(NumSimulParallelotope,1)-.5);
U_2 = Mu(2)+2*s(2)*(rand(NumSimulParallelotope,1)-.5);
U=[U_1 U_2];

% keep ony observations within the ellipsoid
Zsquare = sum((U*inv(Sigma)).*U,2);
Keep=find(Zsquare<=1);
X=U(Keep,:);
NumSimulations=size(X,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

figure 
% marginals
NumBins=round(10*log(NumSimulations));

subplot('Position',[.05 .3 .2 .6]) 
[n,D]=hist(X(:,2),NumBins);
barh(D,n,1);
[y_lim]=get(gca,'ylim')
set(gca,'xtick',[])
grid on

subplot('Position',[.3 .05 .6 .2]) 
[n,D]=hist(X(:,1),NumBins);
bar(D,n,1);
[x_lim]=get(gca,'xlim')
set(gca,'ytick',[])
grid on

% scatter plot
subplot('Position',[.3 .3 .6 .6]) 
h=plot(X(:,1),X(:,2),'.');
set(gca,'xlim',x_lim,'ylim',y_lim)
grid on

% histogram (~rescaled pdf)
NumBins3d=round(sqrt(NumSimulations)/5);
figure
hist3(X,[NumBins3d NumBins3d]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute summary statistics (analytical and sample-based)
Expected_Value=Mu'
Covariance=1/4*Sigma;
Standard_Deviation=sqrt(diag(Covariance))'
Correlation=Covariance(1,2)/(sqrt(Covariance(1,1)*Covariance(2,2)))

Sample_Mean=mean(X)
Sample_Covariance=cov(X);
Sample_Standard_Deviation=sqrt(diag(Sample_Covariance))'
Sample_Correlation=Sample_Covariance(1,2)/(sqrt(Sample_Covariance(1,1)*Sample_Covariance(2,2)))