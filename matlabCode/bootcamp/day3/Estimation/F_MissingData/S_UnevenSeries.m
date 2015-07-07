% this script estimates the parameters of a multivariate normal distribution 
% from an unbalanced panel of time series of different length
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=100;   % maPanel length of panel
Panel(1).Index=[1 4];       % series of equal lenght
Panel(2).Index=[5];       
Panel(3).Index=[3 2];     

Panel(1).s=1;       % first observations in each panel of equal length
Panel(2).s=31;       
Panel(3).s=51;       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build uneven panel
J=length(Panel);
N=0;
for j=1:J
    N_j=length(Panel(j).Index);
    Panel(j).N=N_j;
    N=N+N_j;
end

E=rand(N,1)
sdevs=rand(N,1);
r=.5; % correlation
Corr=(1-r)*eye(N)+r*ones(N,N);
V=diag(sdevs)*Corr*diag(sdevs)

Z=mvnrnd(E,V,T); % generate multivariate normal time series
R=NaN(T,N);
for j=1:J
    R(Panel(j).s:end,Panel(j).Index)=Z(Panel(j).s:end,Panel(j).Index);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[E_hat,V_hat]=UnevenSeriesEstimator(R)