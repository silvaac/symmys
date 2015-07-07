function [mu,sigma_square]=LognMomentsToParameters(e,v)

% this function inverts the formulas (1.98)-(1.99) in the book 
% "Risk and Asset Allocation" - Springer.

sigma_square=log(1+v/(e^2));
mu=log(e)-sigma_square/2;
