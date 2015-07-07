% this script represents covariance and expected value in terms 
% of the location-dispersion ellipsoid in a lognormal example
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters

Mu=[0.15 0.2]';
s=[0.25 0.25];
r=-.9;

NumSimulations=20000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=[s(1)^2     r*s(1)*s(2)
    r*s(1)*s(2)    s(2)^2];

% generate sample
Y = mvnrnd(Mu,Sigma,NumSimulations);
Z = exp(Y);

m=mean(Z)';
s=std(Z)';

% shifted log-normal (covariance = correlation)
X=(Z-ones(NumSimulations,1)*m')./(ones(NumSimulations,1)*s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute location and dispersion (analytical)
Expected_Value=(exp(Mu+(1/2)*diag(Sigma))-m)./s;
Covariance=diag(1./s)*exp(Mu+(1/2)*diag(Sigma))*exp(Mu+(1/2)*diag(Sigma))'.*(exp(Sigma)-1)*diag(1./s);
Standard_Deviation=sqrt(diag(Covariance))
Correlation=Covariance(1,2)/(Standard_Deviation(1)*Standard_Deviation(2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure 
NumBins=round(10*log(NumSimulations));

% scatter plot
subplot('Position',[.3 .3 .6 .6]) 
h=plot(X(:,1),X(:,2),'.');
hold on 
Scale=2;
PlotEigVectors=1;
PlotSquare=1;
TwoDimEllipsoid(Expected_Value,Covariance,Scale,PlotEigVectors,PlotSquare)
[x_lim]=get(gca,'xlim');
[y_lim]=get(gca,'ylim');
grid on

% marginals
subplot('Position',[.05 .3 .2 .6]) 
[n,D]=hist(X(:,2),NumBins);
barh(D,n,1);
set(gca,'ylim',y_lim)
set(gca,'xtick',[])
grid on

subplot('Position',[.3 .05 .6 .2]) 
[n,D]=hist(X(:,1),NumBins);
bar(D,n,1);
set(gca,'xlim',x_lim)
set(gca,'ytick',[])
grid on

