function [ExpectedValue,Volatility, Composition] = EfficientFrontier(NumPortf, Estimate, Constr)
% This function returns the NumPortf x 1 vector expected returns,
%                       the NumPortf x 1 vector of volatilities and
%                       the NumPortf x NumAssets matrix of compositions
% of NumPortf efficient portfolios whose expected returns are equally spaced along the whole range of the efficient frontier
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

warning off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine return of minimum-risk portfolio
FirstDegree=0*Estimate.Exp;
SecondDegree=Estimate.Cov;
Aeq=Constr.Aeq;
beq=Constr.beq;
A=Constr.Aleq;          % no-short constraint
b=Constr.bleq;          % no-short constraint
NumAssets=size(SecondDegree,2);       
x0=1/NumAssets*ones(NumAssets,1);
MinVol_Weights = quadprog(SecondDegree,FirstDegree,A,b,Aeq,beq,[],[],x0);
MinVol_Return=MinVol_Weights'*Estimate.Exp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine return of maximum-return portfolio
[MaxRet_Return,MaxRet_Index]=max(Estimate.Exp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slice efficient frontier in NumPortf equally thick horizontal sectors in the upper branch only
Step=(MaxRet_Return-MinVol_Return)/(NumPortf-1);
TargetReturns=[MinVol_Return : Step : MaxRet_Return];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the NumPortf compositions and risk-return coordinates of the optimal allocations relative to each slice

% start with min vol portfolio
Composition=MinVol_Weights';
Volatility=sqrt(MinVol_Weights'*Estimate.Cov*MinVol_Weights);
ExpectedValue=MinVol_Weights'*Estimate.Exp;

for i=2:NumPortf-1

    % determine least risky portfolio for given expected return
    AEq=[ones(1,NumAssets);
        Estimate.Exp'];
    bEq=[1
        TargetReturns(i)];
    Weights = quadprog(SecondDegree,FirstDegree,A,b,AEq,bEq,[],[],x0)';
    Composition=[Composition
        Weights];
    Volatility=[Volatility
        sqrt(Weights*Estimate.Cov*Weights')];
    ExpectedValue=[ExpectedValue
        Weights*Estimate.Exp];
end
% add max ret portfolio
Weights=zeros(1,NumAssets);
Weights(MaxRet_Index)=1;
Composition=[Composition
    Weights];
Volatility=[Volatility
    sqrt(Weights*Estimate.Cov*Weights')];
ExpectedValue=[ExpectedValue
    Weights*Estimate.Exp];
