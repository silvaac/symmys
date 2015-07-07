% this script shows that lower-order dominance implies higher-order dominance
% compared are two normal distributions
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear;  close all;  clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% order of dominance
Q=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% introduce two discretized pdfs to be compared, normal in this case
mu_A=1;
sig_A=.1;

mu_B=0;
sig_B=.1;

% set up grid for x-axis evaluations of discretized pdfs 
N=2^16;
J=10^6;
dd=normrnd(mu_A,sig_A,1,J);
uu=normrnd(mu_B,sig_B,1,J);
Hi=max([dd uu]);
Lo=min([dd uu]);
h=(Hi-Lo)/(N-1);
X=[Lo+h : h : Hi]';

% discretized pdfs 
Iq_A = 1/h*(normcdf(X+h/2,mu_A,sig_A)-normcdf(X-h/2,mu_A,sig_A));
Iq_B = 1/h*(normcdf(X+h/2,mu_B,sig_B)-normcdf(X-h/2,mu_B,sig_B));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check dominance
for q=2:Q+1
    Iq_A=h*cumsum(Iq_A);
    Iq_B=h*cumsum(Iq_B);
end

Result=['No dominance up to order ' num2str(Q)];
Condition_AdomB=prod(0+(Iq_A<=Iq_B));
if Condition_AdomB
    Result=['A order-' num2str(Q) ' dominates B' ];
end
Condition_BdomA=prod(0+(Iq_B<=Iq_A));
if Condition_BdomA
    Result=['B order-' num2str(Q) ' dominates A' ];
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure
plot(X,Iq_A)
hold on
plot(X,Iq_B,'r')
grid on
legend('A','B')
title(Result)