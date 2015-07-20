clc; clear; close all;

% inputs
tau_tilde=1/52;  % estimation period expressed in years
sig=.4;
P_T=1;

Nsim=10^5;

tau=1/252; tauName = '1 day';           % times to horizon expressed in years
% tau=1/52;  tauName = '1 week';
%tau= 1/12;  tauName = '1 month';
% tau=1; tauName = '1 year';
%  tau=2;   tauName = '2 years'; 


% exact simulation of horizon prices
C_tau=normrnd(0,sig*sqrt(tau),Nsim,1);
P_Ttau=P_T*exp(C_tau);

% compute analytical pdf
p_lo=min(P_Ttau);
p_hi=max(P_Ttau);

p=p_lo : (p_hi-p_lo)/1000 : p_hi; 
%p=linspace(p_lo,p_hi,1000);

m=log(P_T);
s=sig*sqrt(tau);
f=lognpdf(p,m,s);


% compute approximate pdf
f_approx=normpdf(p,P_T,P_T*sig*sqrt(tau));

% plots
figure
Nbins=round(10*log(Nsim));
[f_hist,p_hist]=hist(P_Ttau,Nbins);
f_hist=f_hist/sum((p_hist(2)-p_hist(1))*f_hist); %normalize area under the hist to 1, to make it comparable with a pdf
bar(p_hist,f_hist,'Facecolor',[.7 .7 .7],'Edgecolor',[.5 .5 .5])

hold on
plot(p,f,'r','linewidth',1.5)

hold on
plot(p,f_approx,'b','linewidth',1.5)
xlabel('Price at the horizon')
ylabel('pdf')
title(['Time to horizon \tau = ' tauName])
legend('full Monte Carlo','analytical','Taylor approx')
