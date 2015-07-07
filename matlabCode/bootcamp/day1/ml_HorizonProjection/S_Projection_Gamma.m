% this script performs the horizon projection of a gamma-distributed invariant
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% a detailed proof of the steps below can be found 
% at www.symmys.com > Book > Dowloads > Technical Appendices > Ch 3.7 "Numerical Market Projection"

clc; close all; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
csi=0;
nu=4;
sig=1;

% horizon
T=10;

% coarseness level
N=2^16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up grid
a=-norminv(10^(-15),0,sqrt(T));
h=2*a/N;
Xi=[-a+h : h : a]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=csi+nu*sig*sig;
s=sqrt(2*nu)*sig*sig;
M=(csi-m)/s;

% discretized initial pdf (standardized)
aa=nu/2;  % conversion to matlab notation for inputs
bb=2*sig*sig/s;

f = 1/h*(gamcdf(Xi+h/2-M,aa,bb)-gamcdf(Xi-h/2-M,aa,bb));
f(N) = 1/h*(gamcdf(-a+h/2-M,aa,bb)-gamcdf(-a-M,aa,bb) +...
    gamcdf(a-M,aa,bb)-gamcdf(a-h/2-M,aa,bb));

% discretized characteristic function
Phi=fft(f);                     

% projection of discretized characteristic function
r=[0:(N-1)]';
Signs=(-1).^(r*(T-1));   
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