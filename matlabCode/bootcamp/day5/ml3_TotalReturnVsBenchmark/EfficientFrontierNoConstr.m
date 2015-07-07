function [ExpectedValue,Variance, Composition] = EfficientFrontierNoConstr(Target_ExpectedValues , Covariance, ExpectedValues,Current_Prices,Budget)

warning off;
NumAssets=size(Covariance,2);

% determine exp value of minimum-variance portfolio
FirstDegree=zeros(NumAssets,1);
SecondDegree=Covariance;
Aeq=Current_Prices';
beq=Budget;
A=[];
b=[];
x0=Budget/NumAssets*ones(NumAssets,1);
MinVol_Allocation = quadprog(SecondDegree,FirstDegree,A,b,Aeq,beq,[],[],x0);
MinVol_ExpVal=MinVol_Allocation'*ExpectedValues;

Composition=[];
Variance=[];
ExpectedValue=[];
options = optimset('TolX',1e-15,'Tolfun',1e-15);
for i=1:length(Target_ExpectedValues)
  % determine initial condition
  [Min_ExpectedValue, IndexMin]=min(ExpectedValues);
  [Max_ExpectedValue, IndexMax]=max(ExpectedValues);
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
  Allocation = quadprog(SecondDegree,FirstDegree,A,b,AEq,bEq,[],[],Allocation_0,options)';
  Composition=[Composition 
    Allocation];
  Variance=[Variance
    Allocation*Covariance*Allocation'];
  ExpectedValue=[ExpectedValue
    Allocation*ExpectedValues];
end
