% this script evaluates a robust M-estimator (Huber) of location and scatter 
% by computing replicability, loss, error, bias and inefficiency
% over a stress-test set of correlation values under the multivariate t assumption 
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci
% WARNING: set NumSimulations to ~100 for a quick and dirty result, otherwise it might take a couple of hours

clear; close all; clc; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=4;  % number of joint variables
T=52;  % number of observations in time series

Nu=8;                     % true number of degrees of freedom
Mu=zeros(N,1);            % true expected value
sig=ones(N,1);            % true scatter matrix
Min_Theta=0; Max_Theta=.9; Steps=7; % stress-test the overall correlation of the Student t market

NumSimulations=2000;       % test replicability numerically
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stress test replicability
Step=(Max_Theta-Min_Theta)/(Steps-1);
Thetas=[Min_Theta : Step : Max_Theta];

Stress_Loss_M=[]; Stress_Inef2_M=[]; Stress_Bias2_M=[]; Stress_Error2_M=[];
Stress_Loss_S=[]; Stress_Inef2_S=[]; Stress_Bias2_S=[]; Stress_Error2_S=[];
for i=1:Steps % each cycle represents a different stress-test scenario
    CyclesToGo=Steps-i+1

    Theta=Thetas(i);    
    C=(1-Theta)*eye(N)+Theta*ones(N,N);  
    Sigma= diag(sig)*C*diag(sig);
    S=Sigma*Nu/(Nu-2);
    M=Mu;

    M_hats=[]; S_hats=[];
    l=ones(NumSimulations,1);
    u=ones(T,1);
    for n=1:NumSimulations  % each cycle represents a simulation under a given stress-test scenario
        X=u*Mu' + (u*sig').*mvtrnd(C,Nu,T);
        [M_hat,S_hat]=HubertM(X);

        M_hats=[M_hats
            M_hat(1:end)'];
        S_hats=[S_hats
            S_hat(1:end)];
    end

    % loss for M 
    Loss_M = sum(  (M_hats-l*M').^2  ,2);
    % square inefficiency for M
    Inef2_M = std(M_hats,1)*std(M_hats,1)'; 
    % square bias for M
    Bias2_M = sum(  (mean(M_hats)'-M).^2  ); 
    % square error for M
    Error2_M=mean(Loss_M); 

    % loss for S 
    Loss_S = sum(  (S_hats-l*S(1:end)).^2  ,2);
    % square inefficiency for S
    Inef2_S = std(S_hats)*std(S_hats)'; 
    % square bias for S
    Bias2_S = sum( (mean(S_hats)-S(1:end)).^2 ); 
    % square error for S
    Error2_S=mean(Loss_S); 
    
    % store stress test results
    Stress_Loss_M=[Stress_Loss_M Loss_M];
    Stress_Inef2_M=[Stress_Inef2_M Inef2_M];
    Stress_Bias2_M=[Stress_Bias2_M Bias2_M];
    Stress_Error2_M=[Stress_Error2_M Error2_M];
    Stress_Loss_S=[Stress_Loss_S Loss_S];
    Stress_Inef2_S=[Stress_Inef2_S Inef2_S];
    Stress_Bias2_S=[Stress_Bias2_S Bias2_S];
    Stress_Error2_S=[Stress_Error2_S Error2_S];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
save TempDB_StressCorrelation
clear; clc; close all; load TempDB_StressCorrelation
h=PlotEstimatorStressTest(Stress_Loss_M,Stress_Inef2_M,Stress_Bias2_M,...
    Stress_Error2_M,Thetas,'Correlation','M');
h=PlotEstimatorStressTest(Stress_Loss_S,Stress_Inef2_S,Stress_Bias2_S,...
    Stress_Error2_S,Thetas,'Correlation','S');