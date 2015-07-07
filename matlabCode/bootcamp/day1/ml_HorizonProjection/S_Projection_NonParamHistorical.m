% this script performs the horizon projection of an invariant whose
% fully generic distribution is represented in terms of simulations/historical realizations
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% a detailed proof of the steps below can be found 
% at www.symmys.com > Book > Dowloads > Technical Appendices > Ch 3.7 "Numerical Market Projection"

clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=5; % days to the investment horizon 

N=2^15; % coarseness level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input distribution: choose among
% Spot_USD_CHF Spot_USD_EUR Spot_USD_JPY
% USD_SWAP_2Y_rate USD_SWAP_5Y_rate USD_SWAP_10Y_rate

load('Spot_USD_JPY');
% 1-day empirical observations 
X = log(XX(2:end))-log(XX(1:end-1)); 

% T-day empirical observations 
XXX = XX(1:T:end);
X_empT = log(XXX(2:end))-log(XXX(1:end-1)); 

NumSimul=length(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% standardize distribution
m=mean(X);
s=std(X);
Y=(X-m)/s;

% discretized initial pdf (standardized)
a=-norminv(10^(-15),0,sqrt(T));
h=2*a/N;
Xi=[-a+h : h : a]';
x=Xi+h/2;   

NumPerBin=hist(Y,x);
f = 1/h * NumPerBin'/NumSimul;
f(N) = 1/h * (   NumSimul-sum(NumPerBin(1:end-1))   )/NumSimul;

% discretized characteristic function
Phi=fft(f);                     

% projection of discretized characteristic function
Signs=(-1).^([0:N-1]'*(T-1));   
Phi_T=h^(T-1)*Signs.*(Phi.^T);

% horizon discretized pdf (standardized)
f_T=ifft(Phi_T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-standardized initial pdf 
X_Start=m+s*Xi;
f_Start=f/s;
F_Start=h*cumsum(f_Start*s);

% non-standardized horizon pdf and cdf
x_Hor=m*T+s*Xi;
f_Hor=f_T/s;
F_Hor=h*cumsum(f_Hor*s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

subplot(2,1,1)
u=bar(X_Start,f_Start);
grid on
title('starting pdf')

subplot(2,1,2)
u=bar(x_Hor,f_Hor,'r');
hold on
[zz,x_empT]=hist(X_empT,30);
NobsT=length(X_empT);
hT=x_empT(2)-x_empT(1);
f_empT=zz/(hT*NobsT);
%f_empT=zz/(hT*sum(zz));
u=bar(x_empT,f_empT);
grid on 
legend('empirical pdf','projected pdf','location','northwest')


figure

subplot(2,1,1)
u=bar(X_Start,F_Start);
grid on
title('starting cdf')

subplot(2,1,2)
u=bar(x_Hor,F_Hor,'r');
[F_empT,x_empT]=ecdf(X_empT);
hold on
u=plot(x_empT,F_empT);
grid on 
legend('empirical cdf','projected cdf','location','northwest')