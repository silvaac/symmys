% this script implements the Expectation-Maximization (EM) algoritm, which estimates 
% the parameters of a multivariate normal distribution when some observations are randomly missing
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs

N=2;   % number of variables
T=100;  % lenght of time series
DropNum=4; % number of observations dropped

% parameters of multivariate normal time series
M=0*ones(N,1);
sdevs=ones(N,1);
r=0.99;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate multivariate normal time series
C=(1-r)*eye(N)+r*ones(N,N);
S=diag(sdevs)*C*diag(sdevs);
Original_Data= mvnrnd(M,S,T);

% drop data at random
RData=reshape(Original_Data,T*N,1);
DropIndex=unique(ceil(T*N*rand(DropNum,1)));
RData(DropIndex)=NaN;
Data=reshape(RData,T,N);

% run EM
[E_EM, S_EM , Recovered_Data]=EMestimator(Data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display results
Dates=[1:T];
figure
for n=1:N
    Bad_Dates=Dates(isnan(Data(:,n)));
    Good_Dates=Dates(~isnan(Data(:,n)));
    
    subplot(N,1,n)
    h1=plot(Good_Dates,Original_Data(Good_Dates,n),'.');
    hold on
    h2=plot(Bad_Dates,Original_Data(Bad_Dates,n),'.');
    set(h2,'color','g','markersize',15);
    hold on
    h3=plot(Bad_Dates,Recovered_Data(Bad_Dates,n),'.');
    set(h3,'color','r','markersize',14);
    grid on
end
legend([h2 h3],'original data','EM-recovered','location',[0.2 .03 .05 .03])