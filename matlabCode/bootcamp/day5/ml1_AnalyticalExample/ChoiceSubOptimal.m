function Allocation = ChoiceSubOptimal(Market,InvestorProfile)

N=length(Market.CurrentPrices);

% this is the equally weighted allocation ....

% Allocation=1/N*InvestorProfile.Budget./Market.CurrentPrices;

% ...for graphical and didactical purporses we "cheat": we compute instead an allocation that 
% plots a little to the right of the frontier
Perturbation=2;
Between=.4;

ExpectedValues=diag(Market.CurrentPrices)*(1+Market.LinRets_EV);
Covariance=diag(Market.CurrentPrices)*Market.LinRets_Cov*diag(Market.CurrentPrices);
S=inv(Covariance);
A=Market.CurrentPrices'*S*Market.CurrentPrices; 
B=Market.CurrentPrices'*S*ExpectedValues; 
MV_Portf=InvestorProfile.Budget*S*Market.CurrentPrices/A;
MV_ExpectedValue=MV_Portf'*ExpectedValues; 
SR_Portf=InvestorProfile.Budget*S*ExpectedValues/B;
SR_ExpectedValue=SR_Portf'*ExpectedValues; 

Target_ExpectedValue=MV_ExpectedValue+Between*(SR_ExpectedValue-MV_ExpectedValue);
FirstDegree=zeros(N,1);
SecondDegree=(Covariance+Covariance')/2;
A=[];
b=[];
Aeq=[Market.CurrentPrices'
  ExpectedValues'];
beq=[InvestorProfile.Budget
  Target_ExpectedValue];
Allocation_0=InvestorProfile.Budget./(N*Market.CurrentPrices);
options = optimset('TolX',1e-15,'Tolfun',1e-15);
Allocation = quadprog(SecondDegree,FirstDegree,A,b,Aeq,beq,[],[],Allocation_0,options);

Allocation(3:end)=Allocation(3:end)*(1+Perturbation);
M_1 = [Market.CurrentPrices(1:2)'
  ExpectedValues(1:2)'];
M_2 = [Market.CurrentPrices(3:end)'
  ExpectedValues(3:end)'];
g = [InvestorProfile.Budget
  Target_ExpectedValue];

Allocation(1:2) = inv(M_1)*(g-M_2*Allocation(3:end));