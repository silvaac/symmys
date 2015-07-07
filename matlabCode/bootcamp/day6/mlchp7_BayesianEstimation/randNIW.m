function [Mu,Sigma,InvSigma]=randNIW(Mu_0,T_0,Sigma_0,nu_0,J)
% this function generates a multivariate i.i.d. sample of lenght J from the
% normal-inverse-Wishart distribution: Mu|Sigma ~ N(Mu_0,Sigma/T_0)
%                                      inv(Sigma) ~ W(Nu_0,inv(Sigma_0)/Nu_0)
% See Ch.7 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci

VecIndex=[];
N=length(Mu_0);
for n=1:N
    VecIndex=[VecIndex (n-1)*N+[n:N]];
end

Phi=inv(Sigma_0)/nu_0;
Mu=[];
Sigma=[];
InvSigma=[];
for j=1:J
    Inv_Sigma=wishrnd(Phi,nu_0);
    InvSigma=[InvSigma
        Inv_Sigma(VecIndex)];
    
    S=inv(Inv_Sigma);
    Sigma=[Sigma
        S(VecIndex)];
    
    M = mvnrnd(Mu_0,S/T_0);
    Mu = [Mu
        M];
end


