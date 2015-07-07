% Script for exogenous loadings and endogenous factors 
% Analysis of residual
clear; clc; close all

%------------------------------
% Input parameters:
N = 4;    % market size
numSimulations = 100000;
%------------------------------

% Parameters for the normal market
mu = .1+.3*rand(N,1);
sigma = .5*mu; 
dd=randn(N,N);
[dd,Corr]=cov2corr(dd*dd');

Sigma = diag(sigma)*Corr*diag(sigma); 

% generate simulations for X
X = MvnRnd(mu,Sigma,numSimulations);

% generate a random vector beta
beta = ones(N,1) + randn(N,1)*0.1;

% compute factor realization by cross-sectional regression
F = (X*beta)/(beta'*beta);
% compute residual
U = X - F*beta';
% correlation of residuals U among themselves and with factors F
corr([F U])