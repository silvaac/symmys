function f = LogTPDF(x,nu,Mu,Sigma)
% this function computes the pdf of the N-variate log-t distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (2.213)


N=length(x);

z2=(log(x)-Mu)'*inv(Sigma)*(log(x)-Mu);
K=(nu*pi)^(-N/2) * gamma((nu+N)/2)/gamma(nu/2) * (  (det(Sigma))^(-.5) );

f = K*(1+z2/nu)^(-(nu+N)/2)/prod(x);



