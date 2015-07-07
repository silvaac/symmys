% this script illustrates PCA dimension reduction in a simple bivariate case
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci


clc; close all; clear;

Mu=[0.1 0.2];
Sig=[.25 .25];
rho_12=.8;

NumSimulations=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=diag(Sig)*[1 rho_12 ;  rho_12 1]*diag(Sig);

% compute sample: the market X is the first entry, and the two factors F1 and F2 are the last two entries
Y=mvnrnd(Mu,Sigma,NumSimulations);
Simul_X=Y;%exp(Y);

% compute recovered variables
Expected_Value=mean(Simul_X)';
Covariance=cov(Simul_X);
[EigenVectors,EigenValues] = pcacov(Covariance);

E_k=EigenVectors(:,1);
Recovered_X=[];
for t=1:NumSimulations
    X=Simul_X(t,:)';
    X_p=Expected_Value+E_k*E_k'*(X-Expected_Value);
   
    Recovered_X=[Recovered_X
        X_p' ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plot original sample
figure
h=plot(Simul_X(:,1),Simul_X(:,2),'.');
% plot the ellipsoid
hold on
Scale=2;
PlotEigVectors=1;
PlotSquare=0;
TwoDimEllipsoid(Expected_Value,Covariance,Scale,PlotEigVectors,PlotSquare)

% plot recovered variables
figure
h=plot(Recovered_X(:,1),Recovered_X(:,2),'.');
% plot the ellipsoid
hold on
Scale=2;
PlotEigVectors=1;
PlotSquare=0;
TwoDimEllipsoid(Expected_Value,Covariance,Scale,PlotEigVectors,PlotSquare)