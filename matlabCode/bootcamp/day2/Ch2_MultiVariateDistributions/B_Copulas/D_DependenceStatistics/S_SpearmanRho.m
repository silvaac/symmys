% this script computes Spearman's rho in simulation
% for a log-normal distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

close all; clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
NumSimul=50000;

% input for bivariate normal distribution
r=.999;
Mu=[1 3]';      % NOTE: this input plays no role in the final output
s=[2 5]';    % NOTE: this input plays no role in the final output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate distribution
Sigma = diag(s)*[1 r;r 1]*diag(s);

Y = mvnrnd(Mu,Sigma,NumSimul);      
X = exp(Y);                         % bi-variate log-normal simulation
U_1=logncdf(X(:,1),Mu(1),s(1)); % grade 1 simulation
U_2=logncdf(X(:,2),Mu(2),s(2)); % grade 2 simulation
U=[U_1 U_2];                        % copula

% generate independent copy
YY = mvnrnd(Mu,Sigma,NumSimul);      
XX = exp(YY);                             % bi-variate log-normal simulation
UU_1=logncdf(XX(:,1),Mu(1),s(1)); % grade 1 simulation
UU_2=logncdf(XX(:,2),Mu(2),s(2)); % grade 2 simulation
UU=[UU_1 UU_2];                        % copula

% generate another independent copy
YYY = mvnrnd(Mu,Sigma,NumSimul);     
XXX = exp(YYY);                             % bi-variate log-normal simulation
UUU_1=logncdf(XXX(:,1),Mu(1),s(1)); % grade 1 simulation
UUU_2=logncdf(XXX(:,2),Mu(2),s(2)); % grade 2 simulation
UUU=[UUU_1 UUU_2];                        % copula

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute rho

% sample-based equivalent representation
Concord_Up=(X(:,1)-XX(:,1)).*(X(:,2)-XXX(:,2))>0;
Concord_Down=(X(:,1)-XX(:,1)).*(X(:,2)-XXX(:,2))<0;
Rho_Sample=3*sum(Concord_Up-Concord_Down)/NumSimul

% sample-based from grades
C=corr(U);
Rho_Grades=C(1,2)


