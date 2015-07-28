% This script illustrates the Entropy Pooling approach,
% see A. Meucci, "Entropy Pooling: Fully Flexible Views and
% Stress-Testing," Risk Magazine, October 2008. (Available at
% www.symmys.com > Research.)
clear; clc; close all

% market simulations
nSimulations = 100000;
B = (rand(nSimulations,1) < .5);
X = B.*normrnd(-1,1,nSimulations,1) + (1-B).*normrnd(1,1,nSimulations,1);

% specify view E{X} = 0.5 and constraint 1'*p = 1.
p_prior = ones(nSimulations,1)/nSimulations;
Aeq = [X'; ones(1,nSimulations)];
beq = [0.5; 1];

% posterior market distribution using the Entropy Pooling approach
p_post = EntropyProg(p_prior,[],[],Aeq,beq);

PlotDistributions(X,p_prior,p_post);
fprintf('prior sample mean = %f\n', mean(X));
fprintf('posterior sample mean = %f\n', X'*p_post);
