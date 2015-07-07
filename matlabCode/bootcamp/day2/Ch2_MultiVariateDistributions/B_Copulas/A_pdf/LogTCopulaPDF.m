function F_U = LogTCopulaPDF(u,nu,Mu,Sigma)
% this function computes the pdf of the copula of the log-t distribution
% at the generic point u in the unit hypercube
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci



N=length(u);
s=sqrt(diag(Sigma));

x=exp(Mu+s.*tinv(u,nu));

z2=(log(x)-Mu)'*inv(Sigma)*(log(x)-Mu);
K=(nu*pi)^(-N/2) * gamma((nu+N)/2)/gamma(nu/2) * (  (det(Sigma))^(-.5) );
Numerator = K*(1+z2/nu)^(-(nu+N)/2)/prod(x);

fs=(1./(x.*s)).*tpdf((log(x)-Mu)./s,nu);
Denominator = prod(fs);

F_U = Numerator/Denominator;


