function L=Lklhd(Theta,X,T,N,VecIndex)
% this function computes the likelihood of an N-variate normal sample of length T
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

M = Theta(1:N);
InvS = zeros(N,N);
InvS(VecIndex) = Theta(N+1:end);
InvS=InvS+InvS'-diag(diag(InvS));

lDSig = log(det(InvS));

LnL=0;
for t=1:T
    x = X(t,:)';
    lnf = -N/2*log(2*pi) + .5*lDSig - .5*(x-M)'*InvS*(x-M);
    LnL = LnL+lnf;
end
L=exp(LnL);


