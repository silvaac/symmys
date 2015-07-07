function F_U = TCopulaPDF(u,nu,Mu,Sigma)

N=length(u);
s=sqrt(diag(Sigma));

x=Mu+s.*tinv(u,nu);

z2=(x-Mu)'*inv(Sigma)*(x-Mu);
K=(nu*pi)^(-N/2) * gamma((nu+N)/2)/gamma(nu/2) * (  (det(Sigma))^(-.5) );
Numerator = K*(1+z2/nu)^(-(nu+N)/2);

fs=tpdf((x-Mu)./s,nu)./s;
Denominator = prod(fs);

F_U = Numerator/Denominator;


