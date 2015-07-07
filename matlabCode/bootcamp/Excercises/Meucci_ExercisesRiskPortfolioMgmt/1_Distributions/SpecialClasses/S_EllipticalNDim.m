% this script decomposes the N-variate normal distribution into its radial and uniform components
% then it uses the uniform component to generate an elliptical
% distribution with location parameter Mu and dispersion parameter Sigma
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=30;
NumSimul=10000;
nu=0.1;
t2=0.04;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mu=rand(N,1);
A=rand(N,N)-.5;
Sigma=A*A';

Y = mvnrnd(zeros(N,1),eye(N),NumSimul);

R = sqrt(sum(Y.*Y,2));     % radial distribution (normal case ~ square root of chi-square with N degrees of freedom)
U = Y./(R*ones(1,N));      % uniform distribution on unit sphere

tau=sqrt(t2);
R_New = lognrnd(nu,tau,NumSimul,1);
X = ones(NumSimul,1)*Mu' + (R_New*ones(1,N)).*(U*A');   % N-variate elliptical distribution 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
% visualize projection on m-n coordinates
m=1;
n=3;
figure
plot(X(:,m),X(:,n),'.')
axis equal
grid on
xlabel(sprintf('X_%d',m))
ylabel(sprintf('X_%d',n))

% visualize n-th marginal
n=4;
NumBins=round(10*log(NumSimul));
figure
hist(X(:,n),NumBins);
grid on
xlabel(sprintf('X_%d',n))
