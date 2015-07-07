clear; clc; close all

% project summary statistics to arbitrary horizons 
% see Meucci, A. (2010)
% "Annualization and General Projection of Skewness, Kurtosis and All Summary Statistics"
% GARP Risk Professional - "The Quant Classroom", pp. 52-54 
% Available as "Quant Nugget 4" at http://ssrn.com/abstract=1635484

% last version of this code available at MATLAB Central File Exchange
% http://www.mathworks.com/matlabcentral/fileexchange/25010


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=6;  % focus on first N standardized summary statistics
K=100; % projection horizon

% generate arbitrary distribution
J=100000;  % number of scenarios

Z=randn(J,1); 
X=sin(Z)+log(cos(Z)+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute single-period standardized statistics and central moments
[ga,mu]=SummStats(X,N);

% compute single-period non-central moments
mu_=Central2Raw(mu);

% compute single-period cumulants
ka=Raw2Cumul(mu_);

% compute multi-period cumulants
Ka=K*ka;

% compute multi-period non-central moments
Mu_=Cumul2Raw(Ka);

% compute multi-period central moments
Mu=Raw2Central(Mu_);

% compute multi-period standardized statistics
Ga=Mu;
Ga(2)=sqrt(Mu(2));
for n=3:N
    Ga(n)=Mu(n)/(Ga(2)^n);
end