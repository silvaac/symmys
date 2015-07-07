%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean up workspace
clear; close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define input parameters
NumSimulations=10000; % scalar
sigma_square=6;
ExpX=2;
VarX=7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate t sample with above parameters
mu=ExpX;
nu=2/(1-sigma_square/VarX);
sigma=sqrt(sigma_square);

% built-in generator
X_a=mu+sigma*trnd(nu,NumSimulations,1);  

% stochastic representation
Y=normrnd(0,sigma,NumSimulations,1);
Z=chi2rnd(nu,NumSimulations,1);
X_b=mu+Y./sqrt(Z/nu);

% grade inversion
U=rand(NumSimulations,1);
X_c=mu+sigma*tinv(U,nu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot histograms
NumBins=round(10*log(NumSimulations));

figure % open new figure

subplot(3,1,1)
hist(X_a,NumBins);
grid on
xlabel('built-in generator')
axisLimits = axis;

subplot(3,1,2)
hist(X_b,NumBins);
grid on
xlabel('stoch. representation')
axisLimits(2,:) = axis;

subplot(3,1,3)
hist(X_c,NumBins);
grid on
xlabel('grade inversion')
axisLimits(3,:) = axis;

axisLimits = [min(axisLimits(:,1)),max(axisLimits(:,2)), ...
    min(axisLimits(:,3)),max(axisLimits(:,4))];
subplot(3,1,1), axis(axisLimits);
subplot(3,1,2), axis(axisLimits);
subplot(3,1,3), axis(axisLimits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare empirical quantiles
u=[.01 : .01 : .99];  % range of quantiles (values between zero and one)
q_a=prctile(X_a,u*100);
q_b=prctile(X_b,u*100);
q_c=prctile(X_c,u*100);
figure % open new figure
plot(u,q_a,'b');
hold on 
plot(u,q_b,'g');
hold on 
plot(u,q_c,'r');
grid on
legend('built-in generator','stoch. representation','grade inversion','location','best')
title('quantile of t distribution')
xlabel('Grade')
ylabel('Quantile')