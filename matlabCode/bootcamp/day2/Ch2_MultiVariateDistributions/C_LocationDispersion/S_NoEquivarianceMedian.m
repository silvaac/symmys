% this script is meant to show that the median is not affine equivariant
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
NumSimulations=200;

% distribution
Mu=[0.2 0.1]';  
s=[0.3 0.4];
r=.5;

% invertible affine transformation
m=[.5 .4]';
B=[-1 1.2;.5 -2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate lognormal sample
Sigma=[s(1)^2     r*s(1)*s(2)
    r*s(1)*s(2)    s(2)^2];

Z = mvnrnd(Mu,Sigma,NumSimulations);
X = exp(Z);
Y = ones(NumSimulations,1)*m'+X*B';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute sample means, sample medians and affine transformations thereof
m_X=mean(X)';
q_X=prctile(X,50)';

m_Y=mean(Y);
q_Y=prctile(Y,50)';

m_Y_hat=m+B*m_X;
q_Y_hat=m+B*q_X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure % expected value - sample mean

h=plot(X(:,1),X(:,2),'.');
hold on
h=plot(m_X(1),m_X(2),'.');
set(h,'color','r')

hold on 
h=plot(Y(:,1),Y(:,2),'.');
set(h,'color','g')
hold on
h=plot(m_Y(1),m_Y(2),'.');
set(h,'color','r')
hold on
h=plot(m_Y_hat(1),m_Y_hat(2),'.');
set(h,'color','k')

grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure % median - sample median

h=plot(X(:,1),X(:,2),'.');
hold on
h=plot(q_X(1),q_X(2),'.');
set(h,'color','r')

hold on 
h=plot(Y(:,1),Y(:,2),'.');
set(h,'color','g')
hold on
h=plot(q_Y(1),q_Y(2),'.');
set(h,'color','r')
hold on
h=plot(q_Y_hat(1),q_Y_hat(2),'.');
set(h,'color','k')

grid on