% this script generates a uniform sample on the interior of an ellipsoid in
% N dimensions: X=U*R, where U is a uniform distribution on unit sphere and
% R is the radial component, see "Risk and Asset Allocation" - Springer (2005), by A. Meucci
% This efficient approach is due to Xiaoyu Wang - Dec 2006

clc; clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=2;
Mu=rand(N,1);
A=rand(N,N)-.5;
Sigma=A*A';

NumSimul=10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=ones(1,N);
ll=ones(NumSimul,1);

% 1. Generate U
Z=randn(NumSimul,N);
norm_Z=sqrt(sum(Z.*Z,2));
U=Z./(norm_Z*l);

% 2. Generate R
% the pdf of R is proportional to r^(N-1) therefore the cdf of R is r^N
% we use quantile function of R sample R from uniform simulations
Y=rand(NumSimul,1);
R=Y.^(1/N);

% 3. Generate X
X=ll*Mu'+U.*(R*l)*A';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scatter-plot the first two components
figure
plot(X(:,1),X(:,2),'.')
grid on
axis equal
