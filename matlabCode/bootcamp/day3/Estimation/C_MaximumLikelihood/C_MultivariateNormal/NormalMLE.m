function [Mu,Sigma]=NormalMLE(X)

Mu=mean(X)';

T = size(X,1);
X = X - ones(T,1) * Mu';
Sigma = X'*X/T;

