% This script familiarizes the users with multivariate Bayesian estimation. 
% A normal time series is generated a normal-inverse-Wishart prior is set. 
% The ensuing normal-inverse-Wishart posterior is computed and analyzed numerically and analytically.
% See Ch.7 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear;  close all;  clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ('db_SwapRatesUS');
X=Data(2:end,2:end)-Data(1:end-1,2:end);

% use the following line to specialize the results to N=1
% X(:,[2:end])=[];

NumSimulations=10000; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% history : X_t ~ N(M,S), t=1,...,T
[T,N]=size(X);
Mu_Hat=mean(X)';
Sigma_Hat=cov(X,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prior : Mu|Sigma ~ N(Mu_0,Sigma/T_0)
%         Omega == inv(Sigma) ~ W(Nu_0,inv(Sigma_0)/Nu_0)

% set parameters
Mu_0=0*ones(N,1);
T_0=52;
s=.07*ones(N,1);
r=0.9;
C=(1-r)*eye(N)+r*ones(N,N);
Sigma_0=diag(s)*C*diag(s);
Nu_0=52;

% generate simulations
[Mu_Prior_Simul,Sigma_Prior_Simul,InvSigma_Prior_Simul,VecIndex]=randNIW(T_0,Nu_0,Mu_0,Sigma_0,NumSimulations);

% plot results
PlotNIWMarginals(Mu_Prior_Simul,InvSigma_Prior_Simul,Mu_0,T_0,Sigma_0,Nu_0,'prior');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% posterior : Mu|Sigma ~ N(Mu_1,Sigma/T_1)
%             Omega == inv(Sigma) ~ W(Nu_1,inv(Sigma_1)/Nu_1)

% set parameters
T_1=T_0+T;
Mu_1=(T_0*Mu_0+T*Mu_Hat)/T_1;
Nu_1=Nu_0+T;
Sigma_1=(Nu_0*Sigma_0+T*Sigma_Hat+(T*T_0/(T+T_0))*(Mu_0-Mu_Hat)*(Mu_0-Mu_Hat)')/Nu_1;

% generate simulations
[Mu_Post_Simul,Sigma_Post_Simul,InvSigma_Post_Simul,VecIndex]=randNIW(T_1,Nu_1,Mu_1,Sigma_1,NumSimulations);

% plot results
PlotNIWMarginals(Mu_Post_Simul,InvSigma_Post_Simul,Mu_1,T_1,Sigma_1,Nu_1,'posterior');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute statistics
Mu_CE_Num=mean(Mu_Post_Simul)
Mu_CE_Anal=Mu_1';
Mu_Hat=Mu_Hat'
Mu_0=Mu_0'

Mu_Scatter_Num=cov(Mu_Post_Simul)
Mu_Scatter_Anal=Nu_1/(Nu_1-2)*Sigma_1/T_1;

%
Sigma_CE_Num=mean(Sigma_Post_Simul)
Sigma_CE_Anal=Sigma_1(VecIndex);
Sigma_Hat=Sigma_Hat(VecIndex)
Sigma_0=Sigma_0(VecIndex)

Sigma_Scatter_Num=cov(Sigma_Post_Simul)

%
InvSigma_CE_Num=mean(InvSigma_Post_Simul);
S=inv(Sigma_1);
InvSigma_CE_Anal=S(VecIndex);