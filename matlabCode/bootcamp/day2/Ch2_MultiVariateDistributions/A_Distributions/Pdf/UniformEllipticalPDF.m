function f = UniformEllipticalPDF(x,Mu,Sigma)
% this function computes the pdf of the N-variate 
% uniform distribution on an ellipsoid
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (2.145)


N=length(x);
z2=(x-Mu)'*inv(Sigma)*(x-Mu);

if z2<=1
   K=(pi)^(-N/2) * gamma((N+2)/2) * (  (det(Sigma))^(-.5) );
   f=K;
else
   f=0;	
end



