clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
N=5; % market dimension
K=2; % factors dimension
J=100000; % numbers of simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define true hidden loadings B 
B=rand(N,K)-.5; 
B=B/sqrt(1.5*max(max(B*B')));

% define true hidden idiosyncratic variances 
D=ones(N,1)-diag(B*B');

% define true hidden global covariance
S=B*B'+diag(D);

% generate normal variables with matching moments
X = MvnRnd(zeros(N,1),S,J);

% recover loadings B_, idiosyncratic variances D_ and factors F_ by factor analysis
[B_,D_,T_,stats,F_] = factoran(X,K);

% factor analysis recovers the structure exactly however...
S_=B_*B_'+diag(D_);
Match=1-max(max(abs((S-S_)./S)))

% ...the systematic+idiosyncratic decomposition is NOT recovered
U_=X-F_*B_'; % compute residuals
S_U=corr(U_) % compute correlations