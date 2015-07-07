function F_U = LognormalCopulaPDF(u,Mu,Sigma)
% this function computes the pdf of the copula of the lognormal distribution
% at the generic point u in the unit hypercube
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci


N=length(u);
s=sqrt(diag(Sigma));

x=logninv(u,Mu,s);

Numerator = (2*pi)^(-N/2) * (  (det(Sigma))^(-.5) ) /prod(x) * exp(-.5*(log(x)-Mu)'*inv(Sigma)*(log(x)-Mu));

fs=lognpdf(x,Mu,s);
Denominator = prod(fs);

F_U = Numerator/Denominator;


