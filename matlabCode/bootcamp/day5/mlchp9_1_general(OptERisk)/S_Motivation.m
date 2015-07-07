% this script highlights the issue of estimation risk: the "optimal" frontier estimated with
% naive estimators changes wildly across different time series realizations
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clear; clc; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% market
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=10;
T=52;

MinVol=.05;
MaxVol=.4;
r=.3;
rp=2.5;

Corr=(1-r)*eye(N)+r*ones(N,N);
Step=(MaxVol-MinVol)/(N-1);
Vol=[MinVol : Step : MaxVol]';
Cov=diag(Vol)*Corr*diag(Vol);
ExpVal=rp*Cov*ones(N,1)/N;

R=mvnrnd(ExpVal,Cov,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% investment constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Constr.Aeq=ones(1,N);
Constr.beq=1;
Constr.Aleq=[eye(N)
    -eye(N)];
Constr.bleq=[ones(N,1)
    zeros(N,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% knowledge of the market
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Knowledge.Cov=Cov;
Knowledge.ExpVal=ExpVal;

[e,v,x]=EfficientFrontier(Knowledge, Constr);
TrueEffFront.e=e;
TrueEffFront.v=v;
TrueEffFront.x=x;
PlotFrontier(TrueEffFront)
title('true frontier','fontweight','bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimation of the market
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Estimate.Cov=cov(R);
Estimate.ExpVal=mean(R)';

[e_hat,v_hat,x_hat]=EfficientFrontier(Estimate, Constr);
EstEffFront.e=e_hat;
EstEffFront.v=v_hat;
EstEffFront.x=x_hat;
PlotFrontier(EstEffFront)
title('naive estimation frontier','fontweight','bold')