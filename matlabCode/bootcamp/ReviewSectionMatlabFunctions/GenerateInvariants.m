function [X,F]=GenerateInvariants(mu_x,sig_x,nu_f,sig_f,nu_w,Sigma_w,J)

diag_W=0;
for s=1:nu_w
    Z=mvnrnd([0;0],Sigma_w,J);
    diag_W=diag_W+Z.*Z;
end
a_w=nu_w/2;
b_w=2*diag(Sigma_w);
a_f=nu_f/2;
b_f=2*sig_f^2;
U_x=gamcdf(diag_W(:,1),a_w,b_w(1));
X=logninv(U_x,mu_x,sig_x);
U_f=gamcdf(diag_W(:,2),a_w,b_w(2));
F=gaminv(U_f,a_f,b_f);