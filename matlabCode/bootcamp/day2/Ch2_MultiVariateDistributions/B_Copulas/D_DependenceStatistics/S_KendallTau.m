% this script simulates and computes analytically Kendall's tau 
% for a log-normal distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci


close all; clc; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
NumSimul=1000;

% input for bivariate normal distribution
r=-.99;
Mu=[1 3]';      % NOTE: this input plays no role in the final output
s=[2 5]';    % NOTE: this input plays no role in the final output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate bivariate normal distribution
Sigma = diag(s)*[1 r;r 1]*diag(s);

Y = mvnrnd(Mu,Sigma,NumSimul);      
X = exp(Y);                         % bi-variate log-normal simulation
U_1=logncdf(X(:,1),Mu(1),s(1)); % grade 1 simulation
U_2=logncdf(X(:,2),Mu(2),s(2)); % grade 2 simulation
U=[U_1 U_2];                        % copula

% generate independent copy
YY = mvnrnd(Mu,Sigma,NumSimul);      % bi-variate normal simulation
XX = exp(YY);                             % bi-variate log-normal simulation
UU_1=logncdf(XX(:,1),Mu(1),s(1)); % grade 1 simulation
UU_2=logncdf(XX(:,2),Mu(2),s(2)); % grade 2 simulation
UU=[UU_1 UU_2];                        % copula

% sample-based equivalent representation
Concord_Up=(X(:,1)-XX(:,1)).*(X(:,2)-XX(:,2))>0;
Concord_Down=(X(:,1)-XX(:,1)).*(X(:,2)-XX(:,2))<0;
Tau_Sample=sum(Concord_Up-Concord_Down)/NumSimul

% sample-based integration
Tau_Integr=4*mean(NormalCopulaCDF(U,Mu,Sigma)-1/4)

% analytical
Tau_Anal=2/pi*asin(r)

