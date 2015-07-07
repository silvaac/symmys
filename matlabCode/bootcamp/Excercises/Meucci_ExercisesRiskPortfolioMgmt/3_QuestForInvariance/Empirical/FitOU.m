function [Mu,Th,Sig]=FitOU(Y,tau)

% this function fits a multivariate OU process at estimation step tau
% dY_t=-Th*(Y_t-Mu)*dt+S*dB_t where 
% Th: matrix whose eigenvalues have positive real part and
% dB_t: vector of Brownian motions
% S: any matrix
% output is Mu, Sig=S*S' and Th

% see A. Meucci (2009) 
% "Review of Statistical Arbitrage, Cointegration, and Multivariate Ornstein-Uhlenbeck"
% available at ssrn.com

% Code by A. Meucci, June 2009
% Most recent version available at www.symmys.com > Teaching > MATLAB


[T,N]=size(Y);

X=Y(2:end,:);
F=[ones(T-1,1) Y(1:end-1,:)];
E_XF=X'*F/T;
E_FF=F'*F/T;
B=E_XF*inv(E_FF);

Th=-logm(B(:,2:end))/tau;
Mu=inv(eye(N)-B(:,2:end))*B(:,1);

U=F*B'-X;
Sig_tau = cov(U);

N=length(Mu);
TsT=kron(Th,eye(N))+kron(eye(N),Th);

VecSig_tau=reshape(Sig_tau ,N^2,1);
VecSig = inv(eye(N^2)-expm(-TsT*tau))*TsT*VecSig_tau;
Sig=reshape(VecSig,N,N);