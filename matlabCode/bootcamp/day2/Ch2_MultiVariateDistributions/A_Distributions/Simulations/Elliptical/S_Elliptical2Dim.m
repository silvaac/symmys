% this script decomposes the bivariate normal distribution into its norm and uniform components
% then it uses the uniform component to generate an elliptical
% distribution with location parameter Mu and dispersion parameter Sigma
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mu=[0 1]';
s=[1 1.2]';
r=.4;

NumSimul=10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=diag(s)*[1 r; r 1]*diag(s);
A = chol(Sigma)';
% alternatively use the following:
% [E,L]=pcacov(Sigma);
% A=(E*diag(sqrt(L)));

Y = mvnrnd([0 0]',eye(2),NumSimul);

Norm_Y = sqrt(sum(Y.*Y,2));  
U = Y./(Norm_Y*[1 1]);       % uniform distribution on circle

R_Normal = sqrt(chi2rnd(2,NumSimul,1));  % radial distribution (normal case)
X_Normal = ones(NumSimul,1)*Mu' + (R_Normal*[1 1]).*(U*A');   % bi-variate normal distribution : X = m + RAU

R_Arbitrary = chi2rnd(50,NumSimul,1);      % radial distribution (arbitrary case)
X_Arbitrary = ones(NumSimul,1)*Mu' + (R_Arbitrary*[1 1]).*(U*A');   % bi-variate distribution: X = m + RAU


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumBins=round(10*log(NumSimul));

% normal distribution
figure 

subplot(2,2,1)
plot(U(:,1),U(:,2),'.');
axis equal
grid on
xlabel('uniform on circle U')

subplot(2,2,2)
hist(R_Normal,NumBins);
grid on
xlabel('radial normal')

subplot(2,2,3)
plot(X_Normal(:,1),X_Normal(:,2),'.');
axis equal
grid on
xlabel('joint normal')

% strange distribution
figure

subplot(2,2,1)
plot(U(:,1),U(:,2),'.');
axis equal
grid on
xlabel('uniform on circle U')

subplot(2,2,2)
hist(R_Arbitrary,NumBins);
grid on
xlabel('radial arbitrary')

subplot(2,2,3) % multivariate "strange" distribution
plot(X_Arbitrary(:,1),X_Arbitrary(:,2),'.');
axis equal
grid on
xlabel('joint arbitrary')