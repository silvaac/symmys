%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean up workspace
clear; close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define input parameters
NumSimulations=100000; % scalar
mu_St=0;
s2_St=.1;
nu_St=8;  % NOTE: see how the final results change if you increase nu (=4 is enough)

mu_LN=.1;
s2_LN=.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate Student sample with above parameters
s_St=sqrt(s2_St);
X=mu_St+s_St*trnd(nu_St,NumSimulations,1);

% generate lognormal sample with above parameters
s_LN=sqrt(s2_LN);
Y=lognrnd(mu_LN,s_LN,NumSimulations,1);

% sum samples
Z=X+Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot sample
figure % open new figure
ObservationTime=[1:NumSimulations];
plot(ObservationTime,Z,'.'); 
grid on
xlabel('Simulations')
ylabel('Z')
title('sample vs observation time')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot histogram
NumBins=round(10*log(NumSimulations));
figure % open new figure
hist(Z,NumBins);
grid on
xlabel('Z')
title('sample histogram')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot empirical cdf
figure % open new figure
[f,z]=ecdf(Z);
plot(z,f);
grid on
xlabel('Z')
title('empirical cdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot empirical quantile
u=[.01 : .01 : .99];  % range of quantiles (values between zero and one)
q=prctile(Z,u*100);
figure % open new figure
plot(u,q);
grid on
xlabel('Grade')
ylabel('Quantile')
title('empirical quantile')