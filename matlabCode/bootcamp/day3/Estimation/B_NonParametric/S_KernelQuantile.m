% this script computes the non-parametric kernel estimator of the quantile (value at risk)
% the sample quantile represents the specific case where the bandwidth bw -> 0
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define quantile level
c=.01;

% parameters for t-distributed time series (could be a different distribution)
T=100;
mu=3;
sigma=2;
nu=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define kernel for quantile estimation
th=ceil(c*T); % percentile index
bw=T/100; % bandwidth
sm=normpdf([1:T]',th,bw); % smoother
sm=sm/sum(sm); % normalize smoother

% evaulate kernel estimator of quantile
NumSimul=10000;
Q_sim=zeros(NumSimul,1);
for j=1:NumSimul
    X = mu + sigma*trnd(nu,T,1);
    Sort_X=sort(X);
    Q_sim(j)=Sort_X'*sm;
end

Q_an=mu+sigma*tinv(c,nu);
Bias=abs(mean(Q_sim)-Q_an)
Ineff=std(Q_sim)
Err=sqrt(Bias^2+Ineff^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
NumBins=round(10*log(NumSimul));
hist(Q_sim,NumBins)
hold on
h=plot(Q_an,0,'.');
set(h,'markersize',20,'color','r')