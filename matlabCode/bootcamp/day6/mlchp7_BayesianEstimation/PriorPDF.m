function f = PriorPDF(Theta,Theta_0,c_0,N,VecIndex)
% This function computes the joint pdf of the normal-inverse-Wishart distribution. For notation and theory
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

M=Theta(1:N);
InvS = zeros(N,N);
InvS(VecIndex) = Theta(N+1:end);
InvS=InvS+InvS'-diag(diag(InvS));

M_0=Theta_0(1:N);
InvS_0 = zeros(N,N);
InvS_0(VecIndex) = Theta_0(N+1:end);
InvS_0=InvS_0+InvS_0'-diag(diag(InvS_0));

T_0 = c_0(1);
Nu_0 = c_0(2);

G=InvS*T_0;
Q = (M-M_0)'*G*(M-M_0);
fMcIS = sqrt(det(G)) * exp(-.5*Q);

Sig = InvS_0/Nu_0;
fIS = det(Sig)^(-Nu_0/2) * det(InvS)^((Nu_0-N-1)/2) * exp( -.5*trace(inv(Sig)*InvS) );

f = fMcIS*fIS;
