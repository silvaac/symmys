% this script performs the horizon projection of an invariant whose
% fully generic distribution is represented in terms of simulations/historical realizations
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% a detailed proof of the steps below can be found 
% at www.symmys.com > Book > Dowloads > Technical Appendices > Ch 3.7 "Numerical Market Projection"

clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=2; % multiple of the estimation period to the investment horizon 
N=2^15; % coarseness level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input distribution
NumSimul=100000;
%X = normrnd(0,1,NumSimul,1);
%X = trnd(6,NumSimul,1);
X = lognrnd(0,.4,NumSimul,1);
%X = betarnd(1.1,1,NumSimul,1);

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
u=bar(x_Hor,f_Hor);
hold on
y=normpdf(x_Hor,m*T,s*sqrt(T));
u=plot(x_Hor,y,'r');
grid on 
legend('horizon pdf','CLT limit','location','northwest')

figure

subplot(2,1,1)
u=bar(X_Start,F_Start);
grid on
title('starting cdf')

subplot(2,1,2)
u=bar(x_Hor,F_Hor);
hold on
y=normcdf(x_Hor,m*T,s*sqrt(T));
u=plot(x_Hor,y,'r');
grid on 
legend('horizon cdf','CLT limit','location','northwest')