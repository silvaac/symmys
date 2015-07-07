function f = LognormalPDF(x,Mu,Sigma)
% this function computes the pdf of the N-variate lognormal distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (2.218)



N=length(Mu);
K=(2*pi)^(-N/2) * (  (det(Sigma))^(-.5) ) ;
z2 = (log(x)-Mu)'*inv(Sigma)*(log(x)-Mu);
f = K/prod(x)*exp(-.5*z2);



