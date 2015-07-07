function f = BivStandardTPDF(x,y,r,nu)

det=1-r^2;
z2=(x.^2+y.^2-2*r*x.*y)/det;
K=1/(nu*pi*sqrt(det))*gamma(nu/2+1)/gamma(nu/2);
f = K*((1+z2/nu).^(-nu/2-1));
