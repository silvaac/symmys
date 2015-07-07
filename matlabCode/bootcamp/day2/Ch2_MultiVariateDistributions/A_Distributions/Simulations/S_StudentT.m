% this script simulates a bivariate t distribution 
% it shows its pdf as proxied by the 3-D histogram
% it shows the distribution of a generic linear combination (portfolio) of the two variables
% it computes the summary statistics of the market
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters

Mu=[0.1 1.1]';
s=[0.3 0.2]';
r=-.6;
Nu=1;

LinCombination=[1 -2]';

NumSimulations=100000;
StochRepresentation=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=[s(1)^2     r*s(1)*s(2)
    r*s(1)*s(2)    s(2)^2];
Ones=ones(NumSimulations,1);
if StochRepresentation
    % generate sample by equivalent stochastic representation
    N=mvnrnd(0*Mu,Sigma,NumSimulations);
    Chi2=chi2rnd(Nu,NumSimulations,1);
    X = Ones*Mu' + N./sqrt(Chi2*[1 1]/Nu);
else
    % generate sample by rescaling the built-in generator
    X = Ones*Mu' + (Ones*s').*mvtrnd([1 r;r 1],Nu,NumSimulations);
end

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
if Nu>1
    Expected_Value=Mu'
end
if Nu>2
    Covariance=Nu/(Nu-2)*Sigma;
    Standard_Deviation=sqrt(diag(Covariance))'
    Correlation=r
end

Sample_Mean=mean(X)
Sample_Covariance=cov(X);
Sample_Standard_Deviation=sqrt(diag(Sample_Covariance))'
Sample_Correlation=Sample_Covariance(1,2)/(sqrt(Sample_Covariance(1,1)*Sample_Covariance(2,2)))
