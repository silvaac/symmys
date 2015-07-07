function X=MvnRnd(M,S,J)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function generates normal simulations whose sample moments match the population moments
% inputs are same as MATLAB function "mvnrnd"
% see A. Meucci - "Simulations with Exact Means and Covariances", Risk Magazine, July 2009
% available at http://ssrn.com/abstract=1415699

% Code by A. Meucci. This version June 2009. 
% Last version available at MATLAB central as "Simulations with Exact Means and Covariances"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=length(M);

% generate antithetic variables (mean = 0)
Y = mvnrnd(zeros(N,1),S,J/2);
Y = [Y;-Y];

% compute sample covariance: NOTE defined as "cov(Y,1)", not as "cov(Y)"
S_=cov(Y,1);

% solve Riccati equation using Schur method
H=[zeros(N,N) -S_
   -S zeros(N,N)];
[U_,T_] = schur(H,'real');
U = ordschur(U_,T_,'lhp');
U_lu=U(1:N,1:N);
U_ld=U(N+1:end,1:N);

B=U_ld/U_lu;

% affine transformation to match mean and covariances
X=Y*B+repmat(M',J,1);