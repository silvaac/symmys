%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean up workspace
clear; close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define input parameters
% IMPORTANT: always define the inputs in the beginning in order to make the
% script as flexible as possible
NumSimulations=10000; % scalar
mu=3;
sigma_square=5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate normal sample with above parameters
sigma=sqrt(sigma_square);
X=normrnd(mu,sigma,NumSimulations,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot sample
figure % open new figure
ObservationTime=[1:NumSimulations];
plot(ObservationTime,X,'.'); 
grid on
title('normal sample vs observation time')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot histogram
NumBins=round(10*log(NumSimulations));
figure % open new figure
hist(X,NumBins);
grid on
title('histogram of normal sample')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot empirical cdf
figure % open new figure
[f,x]=ecdf(X);
plot(x,f);
grid on

% plot exact cdf
F=normcdf(x,mu,sigma);
hold on  % superimpose graph
h=plot(x,F);
set(h,'color','r');

legend('empirical','exact','location','best')
title('cdf of normal distribution')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot empirical quantile
u=[.01 : .01 : .99];  % range of quantiles (values between zero and one)
q=prctile(X,u*100);
figure % open new figure
plot(u,q);
grid on

% plot exact quantile
Q=norminv(u,mu,sigma);
hold on 
h=plot(u,Q);
set(h,'color','r');

legend('empirical','exact','location','best')
title('quantile of normal distribution')
xlabel('Grade')
ylabel('Quantile')
