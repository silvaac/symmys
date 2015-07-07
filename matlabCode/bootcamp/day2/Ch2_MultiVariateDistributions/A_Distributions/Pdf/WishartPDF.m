function f = WishartPDF(W,nu,S)
% this function computes the pdf of the NxN-variate Wishart distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (2.224)


N=size(W,1);
Arg=[nu-N+1:nu];
K=2^(nu*N/2)*pi^(N*(N-1)/4)*prod(gamma(Arg/2));

f = 1/K * (det(S))^(-nu/2) * (det(W))^((nu-N-1)/2) * exp(-1/2*trace(inv(S)*W));


