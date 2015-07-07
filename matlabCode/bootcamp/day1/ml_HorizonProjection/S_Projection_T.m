% this script performs the horizon projection of a t-distributed invariant
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% a detailed proof of the steps below can be found 
% at www.symmys.com > Book > Dowloads > Technical Appendices > Ch 3.7 "Numerical Market Projection"

clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
nu=5; % t-distribution's degree of freedom
mu=1; % t-distribution's location parameter
sig=4; % t-distribution's scatter parameter
T=22; % multiple of the estimation period to the investment horizon 

% coarseness level
N=2^14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up grid
a=-norminv(10^(-15),0,sqrt(T));
h=2*a/N;
Xi=[-a+h : h : a]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretized initial pdf (standardized)
f = 1/h*(tcdf(Xi+h/2,nu)-tcdf(Xi-h/2,nu));
f(N) = 1/h*(tcdf(-a+h/2,nu)-tcdf(-a,nu) + tcdf(a,nu)-tcdf(a-h/2,nu));

% discretized characteristic function
Phi=fft(f);                     

% projection of discretized characteristic function
Signs=(-1).^([0:N-1]'*(T-1));   
Phi_T=h^(T-1)*Signs.*(Phi.^T);

% horizon discretized pdf (standardized)
f_T=ifft(Phi_T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-standardized initial pdf 
X_Start=mu+sig*Xi;
f_Start=f/sig;
F_Start=h*cumsum(f_Start*sig);

% non-standardized horizon pdf and cdf
x_Hor=mu*T+sig*Xi;
f_Hor=f_T/sig;
F_Hor=h*cumsum(f_Hor*sig);

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
y=normpdf(x_Hor,mu*T,sig*sqrt(nu/(nu-2))*sqrt(T));
u=plot(x_Hor,y,'r');
grid on 
legend('horizon pdf','CLT limit','location','northwest')

break
figure
subplot(2,1,1)
u=bar(X_Start,F_Start);
grid on
title('starting cdf')

subplot(2,1,2)
u=bar(x_Hor,F_Hor);
hold on
y=normcdf(x_Hor,mu*T,sig*sqrt(nu/(nu-2))*sqrt(T));
u=plot(x_Hor,y,'r');
grid on 
legend('horizon cdf','CLT limit','location','northwest')