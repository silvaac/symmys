clc; close all; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters

theta=.1;   % reversion speed
m=.05;       % long term mean 
sigma=.01;   % volatility
T=10^4;     % number of steps
tau=.01;     % discrete time interval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var=sigma^2/2/theta*(1-exp(-2*theta*tau));
sd=sqrt(var);
eps=normrnd(0,sd,T,1);
x=0;

for t=1:T-1
    x(t+1)=m+exp(-theta*tau)*(x(t)-m)+eps(t);
end

figure
plot(x)
title('AR(1) process vs. time')