% this script shows that lower-order dominance implies higher-order dominance
% compared are the normal and the t distributions
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear;  close all;  clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% order of dominance
Q=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% introduce two discretized pdfs to be compared, normal in this case
mu_A=1;
sig_A=2;

mu_B=0;
sig_B=2;
nu=5;

% set up grid for x-axis evaluations of discretized pdfs 
N=2^16;
J=10^6;
dd=normrnd(mu_A,sig_A,1,J);
ee=mu_B+sig_B*trnd(nu,1,J);
Hi=max([dd ee]);
Lo=min([dd ee]);
h=(Hi-Lo)/(N-1);
X=[Lo+h : h : Hi]';

% discretized normal pdf
Iq_A = 1/h*(normcdf(X+h/2,mu_A,sig_A)-normcdf(X-h/2,mu_A,sig_A));

% discretized t pdf
Iq_B = 1/h*(tcdf((X+h/2-mu_B)/sig_B,nu)-tcdf((X-h/2-mu_B)/sig_B,nu));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check dominance
for q=2:Q+1
    Iq_A=h*cumsum(Iq_A);
    Iq_B=h*cumsum(Iq_B);
end

Result=['No dominance up to order ' num2str(Q)];
Condition_AdomB=prod(0+(Iq_A<=Iq_B));
if Condition_AdomB
    Result=['norm order-' num2str(Q) ' dominates t' ];
end
Condition_BdomA=prod(0+(Iq_B<=Iq_A));
if Condition_BdomA
    Result=['t order-' num2str(Q) ' dominates norm' ];
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure
plot(X,Iq_A)
hold on
plot(X,Iq_B,'r')
grid on
legend('norm','t')
title(Result)