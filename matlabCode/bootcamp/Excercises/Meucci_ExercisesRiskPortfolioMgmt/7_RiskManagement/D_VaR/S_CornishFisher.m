% this script computes the quantile (VaR) 
% - analytically, under the lognormal assumption for the objective
% - using the Cornish-Fisher approximation
% See Sec 5.5 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs 

mu=.05; 
sig=.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process data

E_X=exp(mu+sig^2/2);
Sd_X=exp(mu+sig^2/2)*sqrt(exp(sig^2)-1);
Sk_X=sqrt(exp(sig^2)-1)*(exp(sig^2)+2);

c=[.001:.001:.999];
z=norminv(1-c);

Q_CF = E_X + Sd_X*(  z + Sk_X/6*(z.^2-1)  );
Q_true = logninv(1-c,mu,sig);

x=Q_true;
f=lognpdf(x,mu,sig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure
plot(x,f)
grid on
title('pdf')

figure
plot(c,Q_true,'r')
hold on
plot(c,Q_CF);
grid on
legend('true','Cornish-Fisher approx',2)
title('quantile')