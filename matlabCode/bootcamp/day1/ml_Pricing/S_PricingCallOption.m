% this script shows how the goodness of the delta-gamma-gammadelta approximation 
% for the Black-Scholes call option pricing function
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; close all; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
NumSimul=100000;

% investment horizon (square-root rule)
Horizon=22/252;

% call parameters
K=100;
tau=22/252;

% price parameters
drift=.07;
rf=.03;
vol=.40;
Spot=100;

% set the order in the Taylor expansion
Order=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distribution of the compunded-returns for the horizon
mu=(drift-.5*vol^2)*Horizon;
sigma=sqrt(Horizon)*vol;
C_R=normrnd(mu,sigma,NumSimul,1);

% distribution of prices for the horizon
X=Spot*exp(C_R);


% exact call pricing
C_Exact=blsprice(X, K, rf, tau, vol);

% approximate call pricing
m=mean(X); % Taylor's pivot: can pick an arbitrary value, 
           % but choices other than the mean will require higher orders
D0=blsprice(m, K, rf, tau, vol);
D1=blsdelta(m, K, rf, tau, vol);
D2=blsgamma(m, K, rf, tau,vol);
D3=(blsgamma(m*1.01, K, rf, tau,vol)-blsgamma(m*.99, K, rf, tau,vol))/(m*.02);

switch Order
    case 0
    C_Approx = D0 + 0*(X-m);
    case 1
    C_Approx = D0 + D1*(X-m);
    case 2
    C_Approx = D0 + D1*(X-m) + D2/2*((X-m).^2);
    case 3
    C_Approx = D0 + D1*(X-m) + D2/2*((X-m).^2) + D3/6*((X-m).^3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumBins=round(10*log(NumSimul));

figure

% approximate pricing
subplot('Position',[.25 .3 .15 .6]) 
[n,y]=hist(C_Approx,NumBins);
n=n/(NumSimul*(y(2)-y(1)));
h2=barh(y,n);
set(h2,'facecolor','g','edgecolor','g')
[y_lim]=get(gca,'ylim');
set(gca,'xtick',[])
grid on
title(['order ' num2str(Order) ' approx.'])

% true pricing
subplot('Position',[.05 .3 .15 .6]) 
[n,y]=hist(C_Exact,NumBins);
n=n/(NumSimul*(y(2)-y(1)));
h1=barh(y,n);
set(h1,'facecolor','r','edgecolor','r')
set(gca,'ylim',y_lim)
set(gca,'xtick',[])
grid on
title(['exact'])

% underlying market
subplot('Position',[.45 .05 .45 .2]) 
[n,x]=hist(X,NumBins);
n=n/(NumSimul*(x(2)-x(1)));
h3=bar(x,n);
set(h3,'facecolor','b','edgecolor','b')
set(gca,'ytick',[])
[x_lim]=get(gca,'xlim');
grid on

% pricing functions
subplot('Position',[.45 .3 .45 .6]) 
plot(x,blsprice(x, K, rf, tau, vol),'r')
hold on
y=interp1(X,C_Approx,x,'linear','extrap');
plot(x,y,'g')
set(gca,'xlim',x_lim,'ylim',y_lim)
grid on
legend('exact',['order ' num2str(Order) ' approx.'],'location','northwest')