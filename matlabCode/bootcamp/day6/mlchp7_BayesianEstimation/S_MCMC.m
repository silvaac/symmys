% This script familiarizes the user with the MCMC approach to generate
% the Bayesian posterior distribution of the parameters.
% the MCMC results are checked analytically under the normal-inverse-Wishart specification
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci
% see also the technical appendices at symmys.com>book>downloads

clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1; % dimension of the market. If N=1, the model is normal-inverse-gamma
Burn=100;
NumSimul=2000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate historical information
T=100;
Mu = .5+zeros(N,1);
Sigma = eye(N);
X = mvnrnd(Mu,Sigma,T);

Mu_Hat=mean(X)';
Sigma_Hat=cov(X,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify prior
Mu_0 = 1+zeros(N,1);
T_0=100;
Sigma_0 = 1.3*eye(N);
nu_0=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical posterior
T_1=T_0+T;
Mu_1=(T_0*Mu_0+T*Mu_Hat)/T_1;
nu_1=nu_0+T;
Sigma_1=( nu_0*Sigma_0+T*Sigma_Hat+(Mu_0-Mu_Hat)*(Mu_0-Mu_Hat)'/((1/T)+(1/T_0)) )/nu_1;

[Anal_Mu_Post,aa,Anal_InvSigma_Post] = randNIW(Mu_1,T_1,Sigma_1,nu_1,NumSimul);
Anal_Theta_Post=[Anal_Mu_Post Anal_InvSigma_Post];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCMC posterior

VecIndex=[];
N=length(Mu_0);
for n=1:N
    VecIndex=[VecIndex (n-1)*N+[n:N]];
end

inv_Sigma=inv(Sigma);
Theta = [Mu;
    inv_Sigma(VecIndex)'];
inv_Sigma_0=inv(Sigma_0);
Theta_0 = [Mu_0;
    inv_Sigma_0(VecIndex)'];
inv_Sigma_Hat=inv(Sigma_Hat);
Theta_Hat = [Mu_Hat;
    inv_Sigma_Hat(VecIndex)'];
c_0=[T_0
    nu_0];


Theta_1=Theta_0;
MCMC_Theta_Post=[];
tic
while size(MCMC_Theta_Post,1)<Burn+NumSimul
    Xi = mvnrnd(Theta_1,diag((.1*Theta_Hat).^2),1)';
    u = rand;
    Frac = (Lklhd(Xi,X,T,N,VecIndex)*PriorPDF(Xi,Theta_0,c_0,N,VecIndex)) / ...
        (Lklhd(Theta_1,X,T,N,VecIndex)*PriorPDF(Theta_1,Theta_0,c_0,N,VecIndex));
    if u <= min(Frac ,1)
        Theta_1 = Xi;
        MCMC_Theta_Post = [MCMC_Theta_Post; Theta_1'];
    end
end
MCMC_Theta_Post(1:Burn,:)=[];
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
for k=1:length(Theta)
    figure
    Entry=k;     % pick one entry in the parameter set vector
    MCMC=MCMC_Theta_Post(:,Entry);
    Anal=Anal_Theta_Post(:,Entry);

    % plot MCMC-generated and analytically generated sample of entry
    subplot(2,2,1)
    plot([1:NumSimul],MCMC)
    title('MCMC sample')
    y_lim=get(gca,'ylim');
    subplot(2,2,2)
    plot([1:NumSimul],Anal)
    title('analytical sample')
    set(gca,'ylim',y_lim);

    % plot MCMC-generated and analytically generated histogram of entry
    NumBins=round(10*log(NumSimul));
    subplot(2,2,3)
    hist(MCMC,NumBins);
    x_lim=get(gca,'xlim');
    title('MCMC histogram')
    subplot(2,2,4)
    hist(Anal,NumBins);
    set(gca,'xlim',x_lim);
    title('analytical histogram')

end