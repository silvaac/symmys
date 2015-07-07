% this file generates paths from Merton's jump-diffusion model

% see A. Meucci (2009) 
% "Review of Discrete and Continuous Processes in Finance - Theory and Applications"
% available at ssrn.com

% Code by A. Meucci, April 2009
% Most recent version available at www.symmys.com > Teaching > MATLAB


close all; clc; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
ts=[1/252 : 1/252 : 1]; % grid of time values at which the process is evaluated ("0" will be added, too)
J=3; % number of simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate processes

mu=.00;     % deterministic drift
sig=.20; % Gaussian component


l=3.45; % Poisson process arrival rate
a=0; % drift of log-jump
D=.2; % st.dev of log-jump

X=JumpDiffusionMerton(mu,sig,l,a,D,ts,J);    

figure
plot([0 ts],X');
title('Merton jump-diffusion')
xlabel('time')

