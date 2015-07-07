clc; close all; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters

mu=.05;
a=.45;
b=.5;
s=.03;
T=10^3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=randn(T,1);
s2=s^2;
for t=1:T-1
    s2(t+1)=s^2+a*s2(t)+b*(z(t)^2);
    eps(t+1)=mu+sqrt(s2(t+1))*z(t+1);
end

figure
plot(eps);
title('GARCH(1,1) process vs. time')
