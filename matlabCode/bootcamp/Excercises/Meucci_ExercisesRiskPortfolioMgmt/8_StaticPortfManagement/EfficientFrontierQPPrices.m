function [ExpectedValue,Std_Deviation, Composition] = EfficientFrontierQP(NumPortf, Covariance, ExpectedValues,Current_Prices,Budget)
% this function computes the mean-variance efficient frontier by quadratic programming
% see Sec 6.3 in "Risk and Asset Allocation"-Springer (2005), by A. Meucci

warning off;
NumAssets=size(Covariance,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine exp value of minimum-variance portfolio
FirstDegree=zeros(NumAssets,1);
SecondDegree=Covariance;
Aeq=Current_Prices';
beq=Budget;
A=-eye(NumAssets);
b=zeros(NumAssets,1);
x0=Budget/NumAssets*ones(NumAssets,1);
MinVol_Allocation = quadprog(SecondDegree,FirstDegree,A,b,Aeq,beq,[],[],x0);
MinVol_ExpVal=MinVol_Allocation'*ExpectedValues;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine exp value of maximum-expected value portfolio
Max_ExpVal=Budget*max(ExpectedValues./Current_Prices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slice efficient frontier in NumPortf equally thick horizontal sectors in the upper branch only
Target_ExpectedValues=MinVol_ExpVal + [0:NumPortf]'*(Max_ExpVal-MinVol_ExpVal)/(NumPortf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the NumPortf compositions and risk-return coordinates
Composition=[];
Std_Deviation=[];
ExpectedValue=[];

[Min_ExpectedValue, IndexMin]=min(ExpectedValues);
[Max_ExpectedValue, IndexMax]=max(ExpectedValues);

for i=1:NumPortf

    % determine initial condition
    Matrix=[Min_ExpectedValue Max_ExpectedValue
          Current_Prices(IndexMin) Current_Prices(IndexMax)];
    Allocation_0_MinMax=inv(Matrix)*[Target_ExpectedValues(i); Budget];
    
    Allocation_0=zeros(NumAssets,1);
    Allocation_0(IndexMin)=Allocation_0_MinMax(1);
    Allocation_0(IndexMax)=Allocation_0_MinMax(2);
    
    % determine least risky portfolio for given expected return
    AEq=[Aeq
        ExpectedValues'];
    bEq=[beq
        Target_ExpectedValues(i)];
    Allocation = quadprog(SecondDegree,FirstDegree,A,b,AEq,bEq,[],[],Allocation_0)';
    Composition=[Composition 
                    Allocation];
    Std_Deviation=[Std_Deviation
                sqrt(Allocation*Covariance*Allocation')];
    ExpectedValue=[ExpectedValue
                    Target_ExpectedValues(i)];
end
