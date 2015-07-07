% This script projects the distribution of the compounded returns
% from the estimation interval to the investment horizon 
% Then it computes the distribution of prices at the investment horizon 
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% IMPORTANT: the estimation step "FlexM" was written by O. Ledoit and M. Wolf following 
% Ledoit, O., P. Santa-Clara, and M. Wolf, 2003, Flexible multivariate GARCH
% modeling with an application to international stock markets, Review of Economics and Statistics 85, 735-747.

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs

load DB_Equities.mat
Prices=Prices(:,[4 5]);

J=10000; % numbers of MC scenarios
N=size(Prices,2); % numbers of securities 
T=22; %projection horizon 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimation of daily compounded returns distribution

CR=diff(log(Prices)); % extract risk drivers (compounded returns)

demean=1;
eps=.01;
df=500;
[m,A,B,C,H] = FlexM(CR,demean,eps,df);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection to monthly compouded returns distribution

H_=zeros(J,N,N);
for j=1:J
    H_(j,:,:)=H;
end

X_T=zeros(J,N);
for t=1:T
    for j=1:J % WARNING: this loop is for didactical purposes only. In real applications avoid looping
        % compute new return
        e=randn(N,1);
        H=squeeze(H_(j,:,:));
        X=m+chol(H)*e;
        X_T(j,:)=X_T(j,:)+X';
    
        % update for next cycle
        S=X*X';
        H=C+A.*S+B.*H;
        H_(j,:,:)=H;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pricing into linear returns distribution
R=exp(X_T)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display results
figure 
% marginals
NumBins=round(10*log(J));

subplot('Position',[.05 .3 .2 .6]) 
[n,D]=hist(X_T(:,2),NumBins);
barh(D,n,1);
[y_lim]=get(gca,'ylim');
set(gca,'xtick',[])
grid on

subplot('Position',[.3 .05 .6 .2]) 
[n,D]=hist(X_T(:,1),NumBins);
bar(D,n,1);
[x_lim]=get(gca,'xlim');
set(gca,'ytick',[])
grid on

% scatter plot
subplot('Position',[.3 .3 .6 .6]) 
h=plot(X_T(:,1),X_T(:,2),'.');
set(gca,'xlim',x_lim,'ylim',y_lim)
grid on