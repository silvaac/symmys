% this script computes VaR and contributions to VaR from each exposure
% - semi-analytically (under t assumption)
% - naively-empirical (VaR is percentile, contributions are numerical derivative
% - refined-empirical (VaR from kernel smoothing, contributions as expectations
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs 

N=10;  % number of assets
NumSimul=100000;

% market parameters (student t distribution)
Mu=rand(N,1);
A=rand(N,N)-.5;
Sigma=A*A';
nu=70;

% allocation
a=rand(N,1)-.5;

% VaR confidence
c=.95;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate market scenarios
l=ones(NumSimul,1);
diag_s=diag(sqrt(diag(Sigma)));
C=inv(diag_s)*Sigma*inv(diag_s);
X=mvtrnd(C,nu,NumSimul/2);
X=[X         % symmetrize simulations
    -X];
M = l*Mu' + X*diag_s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% semi-analytical approach
VaR_an=Mu'*a+tinv(1-c,nu)*sqrt(a'*Sigma*a)
DVaR_an=Mu+tinv(1-c,nu)*Sigma*a/sqrt(a'*Sigma*a);
ContrVaR_an=a.*DVaR_an;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% empirical naive approach

% compute objective
Psi = M*a;

% compute naive VaR as percentile
VaR_simul=prctile(Psi,(1-c)*100);

% compute naive gradient as perturbation
DVaR_simul=[];
for n=1:N
    e=zeros(N,1); 
    e(n)=.01;
    Psi_up=M*(a+e);    
    VaR_up=prctile(Psi_up,(1-c)*100);

    DVaR_simul=[DVaR_simul
        (VaR_up-VaR_simul)/e(n)];
end

% compute marginal contributions
ContrVaR_simul=a.*DVaR_simul;

BadError=mean((ContrVaR_simul-ContrVaR_an).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% empirical refined approach

% sort market accoring to objective scenarios
[Sort_Psi,Index]=sort(Psi);
Sort_M=M(Index,:);

% define smoothing  kernel for the empirical VaR computation
th=ceil((1-c)*NumSimul); % percentile index
bw=NumSimul/700; % bandwidth
sm=normpdf([1:NumSimul]',th,bw); % smoother
sm=sm/sum(sm); % normalize smoother

% compute VaR as kernel smoohting
VaR_simul_refined=Sort_Psi'*sm

% compute gradient as expectation
DVaR_simul_refined=0*a;
for n=1:N
    DVaR_simul_refined(n)=sm'*Sort_M(:,n);
end
% compute marginal contributions
ContrVaR_simul_refined=a.*DVaR_simul_refined;

GoodError=mean((ContrVaR_simul_refined-ContrVaR_an).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

figure

subplot('Position',[.2 .7 .75 .22])
DisplayCumumlBars(ContrVaR_an)
xlim([0 N+1])
Ylim=get(gca,'ylim');
title(['analytical contributions: error = 0'])
grid on

subplot('Position',[.05 .7 .1 .22])
bar(VaR_an,'r')
set(gca,'ylim',Ylim);
title(['total'])
grid on

% 
subplot('Position',[.2 .37 .75 .22])
DisplayCumumlBars(ContrVaR_simul)
xlim([0 N+1])
set(gca,'ylim',Ylim);
title(['naive contributions: error = ' num2str(BadError)])
grid on

subplot('Position',[.05 .37 .1 .22])
bar(VaR_simul,'r')
set(gca,'ylim',Ylim);
title(['total'])
grid on

%
subplot('Position',[.2 .05 .75 .22])
DisplayCumumlBars(ContrVaR_simul_refined)
xlim([0 N+1])
set(gca,'ylim',Ylim);
title(['refined contributions: error = ' num2str(GoodError)])
grid on

subplot('Position',[.05 .05 .1 .22])
bar(VaR_simul_refined,'r')
set(gca,'ylim',Ylim);
title(['total'])
grid on
