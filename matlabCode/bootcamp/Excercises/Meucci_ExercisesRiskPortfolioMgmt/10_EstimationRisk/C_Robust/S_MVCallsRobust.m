% this script computes the mean-variance frontier of a set of options
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clc; clear; close all;
% make sure the CVX package is installed and ready to use
if ~exist('cvx_begin','file'),
    error(['Please install the CVX package and make sure that the CVX' ...
        ' root directory is in the MATLAB path before running the script.' ...
        ' CVX is available for download at www.stanford.edu/~boyd/cvx.']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% market
load DB
Stock_0=Stock(end,:);
Vol_0=Vol(end,:);
Strike=Stock_0; % ATM strike
Est=mean(diff(Dates))/252; % estimation interval
Hor=2*Est; % investment horizon

% constraints
N=length(Vol_0);
Constr.Aeq=ones(1,N);  % full-investment constraint
Constr.beq=1;
Constr.Aleq=[eye(N)    % no-short
    -eye(N)];
Constr.bleq=[ones(N,1)
    zeros(N,1)];

% simulations
J=10000; % num simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocation process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% quest for invariance
x_Stock=diff(log(Stock));
x_Vol=diff(log(Vol));

% estimation
M=mean([x_Stock x_Vol])';
S=cov([x_Stock x_Vol]);

% projection
M_Hor=M*Hor/Est;
S_Hor=S*Hor/Est;
X=mvnrnd(M_Hor,S_Hor,J);
X_Stock=X(:,1:N);
X_Vol=X(:,N+1:end);

Stock_Hor=repmat(Stock_0,[J 1]).*exp(X_Stock);
Vol_Hor=repmat(Vol_0,[J 1]).*exp(X_Vol);

% pricing
Call_0=[];
Call_Hor=[];
for n=1:N
    Rate=.04;
    Call_0 = [Call_0 blsprice(Stock_0(n), Strike(n), Rate, Expiry(n), Vol_0(n))];
    Call_Hor = [Call_Hor blsprice(Stock_Hor(:,n), Strike(n), Rate, Expiry(n)-Hor, Vol_Hor(:,n))];
end

% mean-variance
L=Call_Hor./repmat(Call_0,[J 1])-1 ;
Estimate.Cov=cov(L);
Estimate.Exp=mean(L)';
NumPortf=40;
[e,v, w] = EfficientFrontier(NumPortf, Estimate, Constr);
PlotFrontier(w,v);
title('MV efficient frontier')

% robust mean-variance
TargetVols=v;
q=3;T=10;
Estimate.ExpDisp=q*Estimate.Cov/T;
[e_rob,v_rob,w_rob]=RobustEfficientFrontier(TargetVols, Estimate, Constr);
PlotFrontier(w_rob,v_rob);
title('Robust MV efficient frontier')
