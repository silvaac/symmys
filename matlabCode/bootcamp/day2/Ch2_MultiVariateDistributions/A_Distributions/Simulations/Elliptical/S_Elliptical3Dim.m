% this script decomposes the trivariate normal distribution into its norm and uniform components
% then it uses the uniform component to generate an elliptical
% distribution with location parameter Mu and dispersion parameter Sigma
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mu=[1 200 4]'; % location parameter
s=[1 .2  4]';  % dispersions
r=[.3 .2  .5]; % correlations

NumSimul=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = mvnrnd([0 0 0]',eye(3),NumSimul);
Norm_Y = sqrt(sum(Y.*Y,2));  
U = Y./(Norm_Y*[1 1 1]);       % uniform distribution on surface of sphere

Sigma=diag(s)*[1 r(1) r(2); r(1) 1 r(3); r(2) r(3) 1]*diag(s);
A = chol(Sigma)';

R_Normal = sqrt(chi2rnd(3,NumSimul,1));
X_Normal = ones(NumSimul,1)*Mu' + (R_Normal*[1 1 1]).*(U*A');   % tri-variate normal distribution 

R_Arbitrary = chi2rnd(100,NumSimul,1);
X_Arbitrary = ones(NumSimul,1)*Mu' + (R_Arbitrary*[1 1 1]).*(U*A');   % tri-variate elliptical distribution 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumBins=round(10*log(NumSimul));

figure  
% multivariate normal distribution
subplot(2,2,1)
plot3(U(:,1),U(:,2),U(:,3),'.');
axis equal
grid on
xlabel('uniform on sphere U')

subplot(2,2,2)
hist(R_Normal,NumBins);
grid on
xlabel('radial normal')

subplot(2,2,3)
plot3(X_Normal(:,1),X_Normal(:,2),X_Normal(:,3),'.');
grid on
xlabel('joint normal')

figure
% multivariate "strange" distribution
subplot(2,2,1)
plot3(U(:,1),U(:,2),U(:,3),'.');
axis equal
grid on
xlabel('uniform on sphere U')

subplot(2,2,2)
hist(R_Arbitrary,NumBins);
grid on
xlabel('radial strange')

subplot(2,2,3) 
plot3(X_Arbitrary(:,1),X_Arbitrary(:,2),X_Arbitrary(:,3),'.');
grid on
xlabel('joint strange')