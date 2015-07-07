function F_U = NormalCopulaCDF(u,Mu,Sigma)

N=length(u);
s=sqrt(diag(Sigma));
x=norminv(u,Mu,s);
F_U = mvncdf(x,Mu,Sigma);


