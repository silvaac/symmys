% this script illustrates regression dimension reduction in a simple bivariate case
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; close all; clear;

% the invariant X corresponds to the first entry, and the two factors F1 and F2 to the last two entries
Mu=[0.1 0.2];
Sig=[.25 .15];
rho_FX=.8;  

NumSimulations=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=diag(Sig)*[1 rho_FX ; rho_FX 1]*diag(Sig);

% compute sample: the market X is the first entry, and the two factors F1 and F2 are the last two entries
Y=mvnrnd(Mu,Sigma,NumSimulations);
Simul_XF=Y;%exp(Y);

% compute recovered variables
Expected_Value=mean(Simul_XF)';
Covariance=cov(Simul_XF);

ExpVal_X=Expected_Value(1);
Covariance_X=Covariance(1,1);
ExpVal_F=Expected_Value(2);
Covariance_F=Covariance(2,2);
Covariance_XF=Covariance(1,2);

B=Covariance_XF*inv(Covariance_F);

Recovered_XF=[];
for t=1:NumSimulations
    f=Simul_XF(t,2);
    Z=ExpVal_X+B*(f-ExpVal_F);
    Recovered_XF=[Recovered_XF
        Z f ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

figure % original variables
% plot random simulations 
hold on
h=plot(Simul_XF(:,2),Simul_XF(:,1),'.');
% plot the recovered variables
hold on
h=plot(Recovered_XF(:,2),Recovered_XF(:,1),'.');
set(h,'color','r');
grid on

