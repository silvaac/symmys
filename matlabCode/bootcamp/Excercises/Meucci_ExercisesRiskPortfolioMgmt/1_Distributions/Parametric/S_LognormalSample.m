%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean up workspace
clear; close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define input parameters
NumSimulations=10000; % scalar
ExpX=3;
VarX=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate lognormal sample with above parameters
[mu,sigma_square]=LognMomentsToParameters(ExpX,VarX);
sigma=sqrt(sigma_square);
X=lognrnd(mu,sigma,NumSimulations,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot sample
figure % open new figure
ObservationTime=[1:NumSimulations];
plot(ObservationTime,X,'.'); 
grid on
title('lognormal sample vs observation time')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot histogram
NumBins=round(10*log(NumSimulations));
figure % open new figure
hist(X,NumBins);
grid on
title('histogram of lognormal sample')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot empirical cdf
figure % open new figure
[f,x]=ecdf(X);
plot(x,f);
grid on

% plot exact cdf
F=logncdf(x,mu,sigma);
hold on  % superimpose graph
h=plot(x,F);
set(h,'color','r');

legend('empirical','exact','location','best')
title('cdf of lognormal distribution')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot empirical quantile
u=[.01 : .01 : .99];  % range of quantiles (values between zero and one)
q=prctile(X,u*100);
figure % open new figure
plot(u,q);
grid on

% plot exact quantile
Q=logninv(u,mu,sigma);
hold on 
h=plot(u,Q);
set(h,'color','r');

legend('empirical','exact','location','best')
title('quantile of lognormal distribution')
xlabel('Grade')
ylabel('Quantile')