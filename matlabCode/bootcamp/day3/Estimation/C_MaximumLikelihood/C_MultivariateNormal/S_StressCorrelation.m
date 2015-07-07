% this script evaluates the ML estimator of location and scatter under the
% multivariate normal assumption by computing replicability, loss, error, bias and inefficiency
% over a stress-test set of correlation values 
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clear; close all; clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=5;  % number of joint variables
T=20;  % number of observations in time series

Mu=ones(N,1);     %  true location parameter
sig=ones(N,1);  % true dispersions
Min_Theta=0; Max_Theta=.9; Steps=7; % stress-test the overall correlation of the normal market

NumSimulations=2000;   %test replicability numerically 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stress test replicability
Step=(Max_Theta-Min_Theta)/(Steps-1);
Thetas=[Min_Theta : Step : Max_Theta];

Stress_Loss_Mu=[]; Stress_Inef2_Mu=[]; Stress_Bias2_Mu=[]; Stress_Error2_Mu=[];
Stress_Loss_Sigma=[]; Stress_Inef2_Sigma=[]; Stress_Bias2_Sigma=[]; Stress_Error2_Sigma=[];
for i=1:Steps % each cycle represents a different stress-test scenario
    CyclesToGo=Steps-i+1

    Theta=Thetas(i);    
    C=(1-Theta)*eye(N)+Theta*ones(N,N);  
    Sigma= diag(sig)*C*diag(sig);

    Mu_hats=[]; Sigma_hats=[];
    l=ones(NumSimulations,1);
    for n=1:NumSimulations  % each cycle represents a simulation under a given stress-test scenario
        X=mvnrnd(Mu,Sigma,T);
        [Mu_hat,Sigma_hat]=NormalMLE(X);

        Mu_hats=[Mu_hats
            Mu_hat(1:end)'];
        Sigma_hats=[Sigma_hats
            Sigma_hat(1:end)];
    end

    % loss for Mu (numerical)
    Loss_Mu = sum(  (Mu_hats-l*Mu').^2  ,2);
    % square inefficiency for Mu
    Inef2_Mu = 1/T*trace(Sigma); % analytical 
    %Inef2_Mu = std(Mu_hats,1)*std(Mu_hats,1)'; % numerical
    % square bias for Mu
    Bias2_Mu = 0; % analytical 
    %Bias2_Mu = sum(  (mean(Mu_hats)'-Mu).^2  ); % numerical
    % square error for Mu
    Error2_Mu=1/T*trace(Sigma); % analytical 
    %Error2_Mu=mean(Loss_Mu);    % numerical

    % loss for Sigma (numerical)
    Loss_Sigma = sum(  (Sigma_hats-l*Sigma(1:end)).^2  ,2);
    % square inefficiency for Sigma
    Inef2_Sigma = 1/T*(1-1/T)*( trace(Sigma^2)+ trace(Sigma)^2 ); % analytical 
    %Inef2_Sigma = std(Sigma_hats)*std(Sigma_hats)'; % numerical
    % square bias for Sigma
    Bias2_Sigma = trace(Sigma^2)/(T^2);  % analytical
    %Bias2_Sigma = sum( (mean(Sigma_hats)-Sigma(1:end)).^2 ); % numerical
    % square error for Sigma
    Error2_Sigma=1/T*(trace(Sigma*Sigma) + (1-1/T)*(trace(Sigma))^2 ); % analytical
    %Error2_Sigma=mean(Loss_Sigma); % numerical
    
    % store stress test results
    Stress_Loss_Mu=[Stress_Loss_Mu Loss_Mu];
    Stress_Inef2_Mu=[Stress_Inef2_Mu Inef2_Mu];
    Stress_Bias2_Mu=[Stress_Bias2_Mu Bias2_Mu];
    Stress_Error2_Mu=[Stress_Error2_Mu Error2_Mu];
    Stress_Loss_Sigma=[Stress_Loss_Sigma Loss_Sigma];
    Stress_Inef2_Sigma=[Stress_Inef2_Sigma Inef2_Sigma];
    Stress_Bias2_Sigma=[Stress_Bias2_Sigma Bias2_Sigma];
    Stress_Error2_Sigma=[Stress_Error2_Sigma Error2_Sigma];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
h=PlotEstimatorStressTest(Stress_Loss_Mu,Stress_Inef2_Mu,Stress_Bias2_Mu,...
    Stress_Error2_Mu,Thetas,'Correlation','Mu');
h=PlotEstimatorStressTest(Stress_Loss_Sigma,Stress_Inef2_Sigma,Stress_Bias2_Sigma,...
    Stress_Error2_Sigma,Thetas,'Correlation','Sigma');