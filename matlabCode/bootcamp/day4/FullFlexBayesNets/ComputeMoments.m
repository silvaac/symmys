function[m,s,C]=ComputeMoments(X,p)

[J,N]=size(X);

m=X'*p;
Sm=X'*(X.*repmat(p,1,N));
S=Sm-m*m';
[s,C]=cov2corr(S);


