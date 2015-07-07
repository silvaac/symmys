% this script decomposes the N-variate normal distribution into its norm and uniform components
% then it uses the uniform component to generate an elliptical
% distribution with location parameter Mu and dispersion parameter Sigma
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=30;
Mu=rand(N,1);
A=rand(N,N)-.5;
Sigma=A*A';

NumSimul=10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y = mvnrnd(zeros(N,1),eye(N),NumSimul);

Norm_Y = sqrt(sum(Y.*Y,2));     
U = Y./(Norm_Y*ones(1,N));      % uniform distribution on unit sphere

nu=.1;
tau=.2;
R = lognrnd(nu,tau,NumSimul,1).^2; % radial distribution (arbitrary)
X = ones(NumSimul,1)*Mu' + (R*ones(1,N)).*(U*A');   % N-variate elliptical distribution 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
% visualize projection on m-n coordinates
m=1;
n=3;
figure
plot(X(:,m),X(:,n),'.')
grid on

% visualize n-th marginal
n=1;
NumBins=round(10*log(NumSimul));
figure
hist(X(:,n),NumBins);
grid on