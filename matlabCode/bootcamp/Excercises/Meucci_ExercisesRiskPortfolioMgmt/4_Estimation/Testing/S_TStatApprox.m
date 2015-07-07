clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
T=5;
J=30000;

mu_x=.1*(rand-.5);
sig_x=.2*rand;

nu_f=ceil(10*rand);
sig_f=.2*rand;

nu_w=ceil(10*rand);
dd=rand(2,2)-.5;
Sigma_w=dd*dd';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute market features in simulation
[X,F]=GenerateInvariants(mu_x,sig_x,nu_f,sig_f,nu_w,Sigma_w,J);

Mu=mean([X F]);
Sigma=cov([X F]);
mu_X=Mu(1);
mu_F=Mu(2);
sig_X=sqrt(Sigma(1,1));
sig_F=sqrt(Sigma(2,2));
rho=Sigma(1,2)/sqrt(Sigma(1,1)*Sigma(2,2));

Alpha=mu_X-mu_F*rho*sig_X/sig_F;
Beta=rho*sig_X/sig_F;
sig=sig_X*sqrt(1-rho^2);

% randomize time series and compute statistics as random variables
mu_X_hat=zeros(1,J);
sig2_X_hat=zeros(1,J);
Alpha_hat=zeros(1,J);
Beta_hat=zeros(1,J);
sig2_a_hat=zeros(1,J);
sig2_b_hat=zeros(1,J);
sig2_U_hat=zeros(1,J);
for j=1:J

    [X,F]=GenerateInvariants(mu_x,sig_x,nu_f,sig_f,nu_w,Sigma_w,T);
    
    % t-stat for mean
    mu_X_hat(j)=mean(X);
    sig2_X_hat(j)=var(X,1);

    % t-stat for regression
    Sigma_XF=[mean(X) mean(X.*F)];
    Sigma_F(1,1)=1;
    Sigma_F(1,2)=mean(F);
    Sigma_F(2,1)=Sigma_F(1,2);
    Sigma_F(2,2)=mean(F.*F);
    inv_Sigma_F=inv(Sigma_F);
    sig2_a_hat(j)=inv_Sigma_F(1,1);
    sig2_b_hat(j)=inv_Sigma_F(2,2);

    AB_hat=Sigma_XF*inv_Sigma_F;
    Alpha_hat(j)=AB_hat(1);
    Beta_hat(j)=AB_hat(2);        
    
    X_=Alpha_hat(j)+Beta_hat(j)*F;
    U=X-X_;
    sig2_U_hat(j)=var(U,1);
end

t_m=(mu_X_hat-mu_X)./sqrt(sig2_X_hat/(T-1));
t_a=(Alpha_hat-Alpha)./sqrt(sig2_a_hat.*sig2_U_hat/(T-2));
t_b=(Beta_hat-Beta)./sqrt(sig2_b_hat.*sig2_U_hat/(T-2));

% display results
figure
NumBins=round(10*log(J));
subplot(3,1,1)
U_m=tcdf(t_m,T-1);
hist(U_m,NumBins)
subplot(3,1,2)
U_a=tcdf(t_a,T-2);
hist(U_a,NumBins)
subplot(3,1,3)
U_b=tcdf(t_b,T-2);
hist(U_b,NumBins)