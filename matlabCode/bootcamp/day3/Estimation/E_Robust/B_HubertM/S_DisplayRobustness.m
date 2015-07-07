% this script shows in terms of the location-dispersion ellipsoid how,
% unlike the sample estimators of location and scatter, Huber M-estimators of 
% location and scatter are not sensitive to outliers and are thus robust
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clear;  close all;  clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate observations
Mu=[0; 0];
r=.99;
sig=[1 1];
T=50;

Outliers=[-3 2];

Sigma=diag(sig)*[1 r; r 1]*diag(sig);
% "good" observations
Sample = mvnrnd(Mu,Sigma,T);
% add "bad" observation(s)
CorruptSample=[Sample
        Outliers];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M-estimators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute estimates 
[Huber_Mean,Huber_Cov]=HubertM(Sample);
[Huber_MeanCorrupt,Huber_CovCorrupt]=HubertM(CorruptSample);

% plot results
figure 

h=plot(Sample(:,1),Sample(:,2),'.');
set(h,'color','b')
hold on
E=TwoDimEllipsoid(Huber_Mean,Huber_Cov,2,0,0);
set(E,'color','b')

hold on
h=plot(Outliers(:,1),Outliers(:,2),'.');
set(h,'color','r')
hold on
E_Star=TwoDimEllipsoid(Huber_MeanCorrupt,Huber_CovCorrupt,2,0,0);
set(E_Star,'color','r')
title('Huber M-robust estimates')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample estimators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute estimates 
Mean_Sample=mean(Sample)';
Cov_Sample=cov(Sample);

Mean_CorruptSample=mean(CorruptSample)';
Cov_CorruptSample=cov(CorruptSample);

% plot results
figure 

h=plot(Sample(:,1),Sample(:,2),'.');
set(h,'color','b')
hold on
E=TwoDimEllipsoid(Mean_Sample,Cov_Sample,2,0,0);
set(E,'color','b')

hold on
h=plot(Outliers(:,1),Outliers(:,2),'.');
set(h,'color','r')
hold on
E_Star=TwoDimEllipsoid(Mean_CorruptSample,Cov_CorruptSample,2,0,0);
set(E_Star,'color','r')

title('sample estimates')
