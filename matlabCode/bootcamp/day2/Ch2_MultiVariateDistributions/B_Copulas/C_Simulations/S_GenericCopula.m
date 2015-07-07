% this script simulates the copula of a generic distribution. 
% In this example we consider the market-crash copula: low correlations in standard cases,
% high correlations when crashes occur, where the occurrence of crashes can be not-fullly correlated
% see A. Meucci - "Beyond Blak-Litterman in Practice: a Five-Step Recipe to Input Views on Non-Normal Markets"
% Risk Magazine (2006), 19, 9, 114-119

clear;  close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: generate simulations from arbitrary distribution (e.g. correlated extreme events)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=3;
NumSimul=100000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% standard market (low correlations)
Mu_s=zeros(N,1);
Sigma_s=eye(N);
X_s = mvnrnd(Mu_s,Sigma_s,NumSimul);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extreme market
Mu_e=zeros(N,1);
r=.8;
Sigma_e=10*((1-r)*eye(N)+r*ones(N,N));
X_e = mvnrnd(Mu_e,Sigma_e,NumSimul);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% occurrence of extreme events

% probabilities of occurrence of extreme events
p = .1*ones(N,1); 
% correlations of occurrence of extreme events
rho=.99;
C=10*((1-rho)*eye(N)+rho*ones(N,N));

Z = mvnrnd(zeros(N,1),C,NumSimul);
I = normcdf(Z);
B = 0 + (I < ones(NumSimul,1)*p');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% final distribution
X = (1-B).*X_s + B.*X_e;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: extract copula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[W,U]=sort(X);
for n=1:N
    x=U(:,n);
    y=[1:NumSimul];
    xi=[1:NumSimul];
    yi = interp1(x,y,xi);
    Copula(:,n)=yi/(NumSimul+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
% marginals
NumBins=round(10*log(NumSimul));

subplot('Position',[.05 .3 .2 .6]) 
[n,D]=hist(Copula(:,2),NumBins);
barh(D,n,1);
[y_lim]=get(gca,'ylim');
set(gca,'xtick',[])
grid on

subplot('Position',[.3 .05 .6 .2]) 
[n,D]=hist(Copula(:,1),NumBins);
bar(D,n,1);
[x_lim]=get(gca,'xlim');
set(gca,'ytick',[])
grid on

% scatter plot
subplot('Position',[.3 .3 .6 .6]) 
h=plot(Copula(:,1),Copula(:,2),'.');
set(gca,'xlim',x_lim,'ylim',y_lim)
grid on
xlabel('grade 1');
ylabel('grade 2');

% 3-d histogram (~rescaled pdf)
NumBins3d=round(sqrt(NumSimul)/5);
figure
hist3(Copula(:,[1 2]),[NumBins3d NumBins3d]);