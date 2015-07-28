function [Nu,Mu,Sigma]=StudentMLE(x)

% This function computes the maximum-likelihood estimate of the
% degrees of freedom Nu, the location vector Mu and the scatter matrix
% Sigma of T i.i.d. observations from the N-variate t
% distribution, organized in the TxN panel x.
% see Section 4.3.1 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci


Tolerance=.01*mean(prctile(x,75)-prctile(x,25)); % robust scale factor
Nus=[1 2 4 7 12 20]; % significant grid of potential degrees of freedom

LL=[]; % log-likelihood
for s=1:length(Nus)
    Nu=Nus(s);
    [Mu,Sigma] = MleRecursionForT(x,Nu,Tolerance);
    LL=[LL LogLik(x,Nu,Mu,Sigma)];
end
%h=plot(Nus,LL);
[a,Index]=max(LL);
Nu=Nus(Index);
[Mu,Sigma] = MleRecursionForT(x,Nu,Tolerance);

