% this script evaluates the ML estimator of location and scatter under the
% multivariate t assumption by computing replicability, loss, error, bias and inefficiency
% over a stress-test set of correlation values 
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci
% WARNING: set NumSimulations to ~100 for a quick and dirty result, otherwise it might take a couple of hours

clear; close all; clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=4;  % number of joint variables
T=52;  % number of observations in time series

Nu=8;                     % true number of degrees of freedom
Mu=zeros(N,1);            % true location parameter
sig=ones(N,1);  % true dispersions
Min_Theta=0; Max_Theta=.9; Steps=7; % stress-test the overall correlation of the Student t market

NumSimulations=2000;       % test replicability numerically (careful, set this number to ~100 for a quick and dirty result)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stress test replicability
Step=(Max_Theta-Min_Theta)/(Steps-1);
Thetas=[Min_Theta : Step : Max_Theta];

Stress_Loss_Nu=[]; Stress_Inef2_Nu=[]; Stress_Bias2_Nu=[]; Stress_Error2_Nu=[];
Stress_Loss_Mu=[]; Stress_Inef2_Mu=[]; Stress_Bias2_Mu=[]; Stress_Error2_Mu=[];
Stress_Loss_Sigma=[]; Stress_Inef2_Sigma=[]; Stress_Bias2_Sigma=[]; Stress_Error2_Sigma=[];
for i=1:Steps % each cycle represents a different stress-test scenario
    CyclesToGo=Steps-i+1

    Theta=Thetas(i);    
    C=(1-Theta)*eye(N)+Theta*ones(N,N);  
    Sigma= diag(sig)*C*diag(sig);

    Nu_hats=[]; Mu_hats=[]; Sigma_hats=[];
    l=ones(NumSimulations,1);
    u=ones(T,1);
    for n=1:NumSimulations  % each cycle represents a simulation under a given stress-test scenario
        X=u*Mu' + (u*sig').*mvtrnd(C,Nu,T);
        [Nu_hat,Mu_hat,Sigma_hat]=StudentMLE(X);

        Nu_hats=[Nu_hats
            Nu_hat];
        Mu_hats=[Mu_hats
            Mu_hat(1:end)'];
        Sigma_hats=[Sigma_hats
            Sigma_hat(1:end)];
    end

    % loss for Nu 
    Loss_Nu = (Nu_hats-Nu).^2;
    % square inefficiency for Nu
    Inef2_Nu = std(Nu_hats,1)^2; 
    % square bias for Nu
    Bias2_Nu = (mean(Nu_hats)-Nu)^2; 
    % square error for Nu
    Error2_Nu=mean(Loss_Nu); 

    % loss for Mu 
    Loss_Mu = sum(  (Mu_hats-l*Mu').^2  ,2);
    % square inefficiency for Mu
    Inef2_Mu = std(Mu_hats,1)*std(Mu_hats,1)'; 
    % square bias for Mu
    Bias2_Mu = sum(  (mean(Mu_hats)'-Mu).^2  ); 
    % square error for Mu
    Error2_Mu=mean(Loss_Mu); 

    % loss for Sigma 
    Loss_Sigma = sum(  (Sigma_hats-l*Sigma(1:end)).^2  ,2);
    % square inefficiency for Sigma
    Inef2_Sigma = std(Sigma_hats)*std(Sigma_hats)'; 
    % square bias for Sigma
    Bias2_Sigma = sum( (mean(Sigma_hats)-Sigma(1:end)).^2 ); 
    % square error for Sigma
    Error2_Sigma=mean(Loss_Sigma); 
    
    % store stress test results
    Stress_Loss_Nu=[Stress_Loss_Nu Loss_Nu];
    Stress_Inef2_Nu=[Stress_Inef2_Nu Inef2_Nu];
    Stress_Bias2_Nu=[Stress_Bias2_Nu Bias2_Nu];
    Stress_Error2_Nu=[Stress_Error2_Nu Error2_Nu];
    Stress_Loss_Mu=[Stress_Loss_Mu Loss_Mu];
    Stress_Inef2_Mu=[Stress_Inef2_Mu Inef2_Mu];
    Stress_Bias2_Mu=[Stress_Bias2_Mu Bias2_Mu];
    Stress_Error2_Mu=[Stress_Error2_Mu Error2_Mu];
    Stress_Loss_Sigma=[Stress_Loss_Sigma Loss_Sigma];
    Stress_Inef2_Sigma=[Stress_Inef2_Sigma Inef2_Sigma];
    Stress_Bias2_Sigma=[Stress_Bias2_Sigma Bias2_Sigma];
    Stress_Error2_Sigma=[Stress_Error2_Sigma Error2_Sigma];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
h=PlotEstimatorStressTest(Stress_Loss_Nu,Stress_Inef2_Nu,Stress_Bias2_Nu,...
    Stress_Error2_Nu,Thetas,'Correlation','Nu');
h=PlotEstimatorStressTest(Stress_Loss_Mu,Stress_Inef2_Mu,Stress_Bias2_Mu,...
    Stress_Error2_Mu,Thetas,'Correlation','Mu');
h=PlotEstimatorStressTest(Stress_Loss_Sigma,Stress_Inef2_Sigma,Stress_Bias2_Sigma,...
    Stress_Error2_Sigma,Thetas,'Correlation','Sigma');