% this script performs all the computations for the comparison of total-return vs. benchmark-driven allocation
% see Sec. 6.6 of "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=5; % number of Securities
Budget=1000;


% generate prices
CurrentPrices=ones(N,1);

% generate exp values and covariances
Correlations=eye(N)
Max_ExpectedValue=10; Min_ExpectedValue=1; Step=(Max_ExpectedValue-Min_ExpectedValue)/(N-1);
ExpectedValues=[Min_ExpectedValue : Step : Max_ExpectedValue]';
StDeviations=4*ExpectedValues;
Covariance=diag(StDeviations)*Correlations*diag(StDeviations);

% generate quasi-MV-efficient benchmark (peturbation of min variance portfolio)
MV_Portf=Budget*inv(Covariance)*CurrentPrices;
Benchmark=MV_Portf.*randn(N,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Norm_Benchmark=Budget*Benchmark/(CurrentPrices'*Benchmark);
Bench_ExpectedValue=Norm_Benchmark'*ExpectedValues;
Bench_Variance=Norm_Benchmark'*Covariance*Norm_Benchmark;
Bench_EOP=0;
Bench_TE2=0;

S=inv(Covariance);
A=CurrentPrices'*S*CurrentPrices;
B=CurrentPrices'*S*ExpectedValues;
C=ExpectedValues'*S*ExpectedValues;
D=A*C-B*B;
Delta_B = Bench_Variance - A/D*Bench_ExpectedValue^2 + 2*Budget*B/D*Bench_ExpectedValue - Budget^2*C/D;

MV_Portf=Budget*S*CurrentPrices/A;
MV_ExpectedValue=MV_Portf'*ExpectedValues;
MV_Variance=MV_Portf'*Covariance*MV_Portf;
MV_EOP=(MV_Portf-Norm_Benchmark)'*ExpectedValues;
MV_TE2=(MV_Portf-Norm_Benchmark)'*Covariance*(MV_Portf-Norm_Benchmark);

Sh_Portf=Budget*S*ExpectedValues/B;
Sh_ExpectedValue=Sh_Portf'*ExpectedValues;
Sh_Variance=Sh_Portf'*Covariance*Sh_Portf;
Sh_EOP=(Sh_Portf-Norm_Benchmark)'*ExpectedValues;
Sh_TE2=(Sh_Portf-Norm_Benchmark)'*Covariance*(Sh_Portf-Norm_Benchmark);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total return objective
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Step=Budget/10;
E=MV_ExpectedValue+[-5*Budget : Step : 5*Budget];

% method 1: the analtical equation of the parabola
TotRetCurve_Es_1=E;
TotRetCurve_Vs_1 = A/D*E.^2 - 2*Budget*B/D*E + Budget^2*C/D;

TotRetCurve_EOPs_1 = E-Bench_ExpectedValue;
TotRetCurve_TE2s_1 = A/D*TotRetCurve_EOPs_1.^2 + Delta_B;

% method 2: from the analytical allocation curve to its coordinates
TotRetCurve_Es_2=[];
TotRetCurve_Vs_2=[];
TotRetCurve_EOPs_2=[];
TotRetCurve_TE2s_2=[];
for i=1:length(E);
    TotRetCurve_Portf=MV_Portf + (E(i)-MV_ExpectedValue)* (Sh_Portf-MV_Portf)/(Sh_ExpectedValue-MV_ExpectedValue);

    TotRetCurve_ExpectedValue=TotRetCurve_Portf'*ExpectedValues;
    TotRetCurve_Es_2=[TotRetCurve_Es_2 TotRetCurve_ExpectedValue];

    TotRetCurve_Variance=TotRetCurve_Portf'*Covariance*TotRetCurve_Portf;
    TotRetCurve_Vs_2=[TotRetCurve_Vs_2 TotRetCurve_Variance];

    TotRetCurve_EOP=(TotRetCurve_Portf-Norm_Benchmark)'*ExpectedValues;
    TotRetCurve_EOPs_2=[TotRetCurve_EOPs_2 TotRetCurve_EOP];

    TotRetCurve_TE2=(TotRetCurve_Portf-Norm_Benchmark)'*Covariance*(TotRetCurve_Portf-Norm_Benchmark);
    TotRetCurve_TE2s_2=[TotRetCurve_TE2s_2 TotRetCurve_TE2];
end

% method 3: from the quadratic programming solution to its coordinates
[TotRetCurve_Es_3,TotRetCurve_Vs_3, Allocations] = EfficientFrontierNoConstr(E,Covariance,ExpectedValues,CurrentPrices,Budget);
TotRetCurve_EOPs_3=[];
TotRetCurve_TE2s_3=[];
for i=1:size(Allocations,1)
    TotRetCurve_Portf=Allocations(i,:)';
    TotRetCurve_EOP=(TotRetCurve_Portf-Norm_Benchmark)'*ExpectedValues;
    TotRetCurve_EOPs_3=[TotRetCurve_EOPs_3 TotRetCurve_EOP];

    TotRetCurve_TE2=(TotRetCurve_Portf-Norm_Benchmark)'*Covariance*(TotRetCurve_Portf-Norm_Benchmark);
    TotRetCurve_TE2s_3=[TotRetCurve_TE2s_3 TotRetCurve_TE2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% benchmark-relative objective
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EOP=[-5*Budget : Step : 5*Budget];

% method 1: the analtical equation of the parabola
BenchRelCurve_EOPs_1 = EOP;
BenchRelCurve_TE2s_1 = A/D*EOP.^2;

BenchRelCurve_Es_1 = EOP + Bench_ExpectedValue;
BenchRelCurve_Vs_1 = A/D*BenchRelCurve_Es_1.^2 - 2*Budget*B/D*BenchRelCurve_Es_1 + Budget^2*C/D + Delta_B;

% method 2: from the analytical allocation curve to its coordinates
BenchRelCurve_Es_2=[];
BenchRelCurve_Vs_2=[];
BenchRelCurve_EOPs_2=[];
BenchRelCurve_TE2s_2=[];
for i=1:length(EOP);
    BenchRelCurve_Portf = Norm_Benchmark + EOP(i)* (Sh_Portf-MV_Portf)/(Sh_ExpectedValue-MV_ExpectedValue);

    BenchRelCurve_ExpectedValue=BenchRelCurve_Portf'*ExpectedValues;
    BenchRelCurve_Es_2=[BenchRelCurve_Es_2 BenchRelCurve_ExpectedValue];

    BenchRelCurve_Variance=BenchRelCurve_Portf'*Covariance*BenchRelCurve_Portf;
    BenchRelCurve_Vs_2=[BenchRelCurve_Vs_2 BenchRelCurve_Variance];

    BenchRelCurve_EOP=(BenchRelCurve_Portf-Norm_Benchmark)'*ExpectedValues;
    BenchRelCurve_EOPs_2=[BenchRelCurve_EOPs_2 BenchRelCurve_EOP];

    BenchRelCurve_TE2=(BenchRelCurve_Portf-Norm_Benchmark)'*Covariance*(BenchRelCurve_Portf-Norm_Benchmark);
    BenchRelCurve_TE2s_2=[BenchRelCurve_TE2s_2 BenchRelCurve_TE2];
end

% method 3: from the quadratic programming solution to its coordinates
c=Budget*10^(-6);  % relatively speaking practically zero
[BenchRelCurve_EOPs_3,BenchRelCurve_TE2s_3, RelAllocations] = EfficientFrontierNoConstr(EOP,Covariance,ExpectedValues,CurrentPrices,c);
BenchRelCurve_Es_3=[];
BenchRelCurve_Vs_3=[];
for i=1:size(RelAllocations,1)
    BenchRelCurve_Portf=RelAllocations(i,:)'+ Norm_Benchmark;

    BenchRelCurve_E=BenchRelCurve_Portf'*ExpectedValues;
    BenchRelCurve_Es_3=[BenchRelCurve_Es_3 BenchRelCurve_E];

    BenchRelCurve_V=(BenchRelCurve_Portf)'*Covariance*(BenchRelCurve_Portf);
    BenchRelCurve_Vs_3=[BenchRelCurve_Vs_3 BenchRelCurve_V];
end

PlotTotalRetsVsRelativeRets