function F_U = TCopulaPDF(u,nu,Mu,Sigma)
% this function computes the pdf of the copula of the t distribution
% at the generic point u in the unit hypercube
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci


N=length(u);
s=sqrt(diag(Sigma));

x=Mu+s.*tinv(u,nu);

z2=(x-Mu)'*inv(Sigma)*(x-Mu);
K=(nu*pi)^(-N/2) * gamma((nu+N)/2)/gamma(nu/2) * (  (det(Sigma))^(-.5) );
Numerator = K*(1+z2/nu)^(-(nu+N)/2);

fs=tpdf((x-Mu)./s,nu)./s;
Denominator = prod(fs);

F_U = Numerator/Denominator;


