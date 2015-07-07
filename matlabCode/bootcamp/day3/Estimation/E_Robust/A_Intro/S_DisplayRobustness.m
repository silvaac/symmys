% this script shows in terms of the location-dispersion ellipsoid how the sample 
% mean and covariance are sensitive to outliers and are thus not robust
% see Sec. 4.5 in "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clear;  close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate observations
Mu=[0; 0];
r=-.90;
sig=[1 1];
T=50;

Outliers=10*rand(1,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=diag(sig)*[1 r; r 1]*diag(sig);
% "good" observations
Sample = mvnrnd(Mu,Sigma,T);
% add "bad" observation(s)
CorruptSample=[Sample
    Outliers];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute sample estimates
Sample_Mean=mean(Sample)';
Sample_Cov=cov(Sample);
Sample_MeanCorrupt=mean(CorruptSample)';
Sample_CovCorrupt=cov(CorruptSample);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results

figure % data
h=plot(Sample(:,1),Sample(:,2),'.');
set(h,'color','b')
E=TwoDimEllipsoid(Sample_Mean,Sample_Cov,2,0,0);
set(E,'color','b','linewidth',1)
hold on
h=plot(Outliers(:,1),Outliers(:,2),'.');
set(h,'color','r')
E=TwoDimEllipsoid(Sample_MeanCorrupt,Sample_CovCorrupt,2,0,0);
set(E,'color','r','linewidth',1)
title('sample estimates')