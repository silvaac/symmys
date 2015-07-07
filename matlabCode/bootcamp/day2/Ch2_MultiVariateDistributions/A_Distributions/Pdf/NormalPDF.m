function f = NormalPDF(x,Mu,Sigma)
% this function computes the pdf of the N-variate normal distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (2.156)



N=length(x);

K=(2*pi)^(-N/2) * (  (det(Sigma))^(-.5) ); 
z2=(x-Mu)'*inv(Sigma)*(x-Mu);
f = K* exp(-.5*z2);


