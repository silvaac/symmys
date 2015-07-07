function F_U = NormalCopulaPDF(u,Mu,Sigma)

N=length(u);
s=sqrt(diag(Sigma));

x=norminv(u,Mu,s);

Numerator = (2*pi)^(-N/2) * (  (det(Sigma))^(-.5) ) * exp(-.5*(x-Mu)'*inv(Sigma)*(x-Mu));

fs=normpdf(x,Mu,s);
Denominator = prod(fs);

F_U = Numerator/Denominator;


