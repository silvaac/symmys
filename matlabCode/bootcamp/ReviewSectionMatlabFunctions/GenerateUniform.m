function X=GenerateUniform(J,N)
% this function generates a uniform sample on the unit hypersphere
% See "Risk and Asset Allocation" - Springer (2005), by A. Meucci
% script by Xiaoyu Wang - Dec 2006
% we decompose X=U*R, where U is a uniform distribution on unit sphere and
% R is a distribution on (0,1) proportional to r^(Dims-1), i.e. the area of
% surface of radius r

l=ones(1,N);

% 1. Generate U
Z=randn(J,N);
normZ=sqrt(sum(Z.*Z,2));
U=Z./(normZ*l);

% 2. Generate R
% the pdf of R is proportional to r^(N-1) therefore the cdf of R is r^N
% we use quantile function of R sample R from uniform simulations
Y=rand(J,1);
R=Y.^(1/N);

% 3. Generate X
X=U.*(R*l);