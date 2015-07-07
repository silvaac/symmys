%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this scripts generates normal simulations whose sample moments match the population moments
% see A. Meucci - "Simulations with Exact Means and Covariances", Risk Magazine, July 2009
% available at www.ssrn.com

% Code by A. Meucci. This version June 2009. 
% Last version available at MATLAB central as "Simulations with Exact Means and Covariances"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs

N=20; % number of risk factors
J=200; % number of simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate desired population moments:

% ...vector of expected values M...
M=rand(N,1)-.5;
% ...covariance matrix S
A=rand(N,N)-.5;
S=A*A';

% generate sample of size J from multivariate normal N(M,S)
X=mvnrnd(M,S,J); % no match between sample and population moments (built-in) function
%X=MvnRnd(M,S,J); % exact match between sample and population moments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute sample moments and errors
M_=mean(X)';
S_=cov(X,1);

Err_M=max(abs(M-M_))/max(abs(M))
Err_S=max(max(abs(S-S_)))/max(max(abs(S)))