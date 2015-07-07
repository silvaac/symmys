% this script computes the cdf of a bivariate normal distribution at a grid of points
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
m=[0.04  0.05];
s=[.2 .25];    
r=.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input generate grid according to natural disribution's limits
C=[1 r;
    r 1];
S=diag(s)*C*diag(s);

Y=mvnrnd(m,S,50);
Grid1=[min(Y(:,1))  :  (max(Y(:,1))-min(Y(:,1)))/25  : max(Y(:,1))];
Grid2=[min(Y(:,2))  :  (max(Y(:,2))-min(Y(:,2)))/25  : max(Y(:,2))];

mu = [1 -1]; Sigma = [.9 .4; .4 .3];
[X1,X2] = meshgrid(Grid1,Grid2);
X = [X1(:) X2(:)];
f = mvncdf(X, m, S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
surf(X1,X2,reshape(f,length(Grid2),length(Grid1)));
