function F_U = NormalCopulaCDF(u,Mu,Sigma)

N=length(u);
s=sqrt(diag(Sigma));
x1=norminv(u(:,1),Mu(1),s(1));
x2=norminv(u(:,2),Mu(2),s(2));

x=[x1 x2];
F_U = mvncdf(x,Mu',Sigma);


