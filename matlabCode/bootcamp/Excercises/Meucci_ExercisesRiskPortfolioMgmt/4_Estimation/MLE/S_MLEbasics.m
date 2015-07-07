% this script performs ML under a non-standard parametric set of distributions
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clear;  close all;  clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
p=.01;
Theta=[[-.04 : .001: -.01] .02 .03];
load DBTimeSeries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check invariance
T=length(i_T);
IIDAnalysis([1:T],i_T)
axisLimits = axis;
text(axisLimits(1:2)*[-0.1,1.1]', axisLimits(3:4)*[0.1,0.9]', ...
    'Invariance Analysis');


% ML-estimate parameters
Store_LL=[];
for s=1:length(Theta)
    s
    theta=Theta(s);
    % compute log-likelihood
    LL=0;
    for t=1:T
        x_t=i_T(t);
        LL=LL+log(fparam(x_t,theta));
    end
    
    Store_LL=[Store_LL LL];
end

[Max,Max_Index]=max(Store_LL);

theta_ML=Theta(Max_Index)

% compute MLE-implied quantile
Q_ML=qparam(p,theta_ML)

% compute sample quantile
Q_NP=prctile(i_T,p*100)