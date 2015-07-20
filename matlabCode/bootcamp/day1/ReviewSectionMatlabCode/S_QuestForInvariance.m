%% EQUITY

clear; clc;
load DB_Equities

pick=20; %consider the 20th stock

%Prices
P=Prices(end-599:end,pick); 
IIDAnalysis(P)

%Total returns
X=P(2:end)./P(1:end-1);
IIDAnalysis(X)

%Price increments
Y=P(2:end)-P(1:end-1); %or Y=diff(P)
IIDAnalysis(Y)

%Squared total returns
Z=X.^2;
IIDAnalysis(Z)

%W
W=P(3:end)-2*P(2:end-1)+P(1:end-2);
IIDAnalysis(W)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIXED INCOME

clear;close all; clc;
load DB_FixedIncome

t2m=ycMaturityYrs;
yields=ycYieldPercent/100;
clear ycMaturityYrs ycYieldPercent

line(t2m,yields(1,:))
line(t2m,yields(2,:),'color','r')
line(t2m,yields(end,:),'color','g')
xlabel('time to maturity')
ylabel('yield to maturity')

%changes in the yield for a specific time to mat
Y=yields(:,t2m==5);
dY=diff(Y);
IIDAnalysis(dY)

%changes in the log-yield for a specific time to mat
X=log(Y);
dX=diff(X);
IIDAnalysis(dX)

perm_idx=randperm(length(dX));
dX1=dX(perm_idx);
 
sample1=dX1(1:length(dX1)/2);
sample2=dX1(length(dX1)/2+1:end);

subplot(1,2,1)
hist(sample1)
subplot(1,2,2)
hist(sample2)

figure()
[f,x]=ecdf(sample1);
plot(x,f)
[f,x]=ecdf(sample2);
hold on
plot(x,f,'r')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DERIVATIVES

clear;close all; clc;
load DB_Derivatives

colormap(bone)
surf(moneyness,days2Maturity,squeeze(impVol(1,:,:)))
xlabel('moneyness (K/S)')
ylabel('time to mat (days)')
zlabel('implied vol.')

pick_t2m=1;
pick_mon=find(moneyness==1);


%weekly changes in implied volatility
X_daily=impVol(:,pick_t2m,pick_mon); %daily time series of imp. vol. for fixed time to mat and moneyness
X_weekly=X_daily(1:5:end);%weekly time series of imp. vol. for fixed time to mat and moneyness
dX_weekly=diff(X_weekly);
IIDAnalysis(dX_weekly)

%weekly changes in log-implied volatility
Y=log(X_weekly);
dY=diff(Y);
IIDAnalysis(dY)


%define variable Z
[T,Mat,Mon]=size(impVol(1:5:end,:,:));
Z=reshape(log(impVol(1:5:end,:,:)),T,Mat*Mon);

% VAR(1) model estimation 
X=Z(2:end,:)'; %dependent variables time series (N x T)
F=Z(1:end-1,:)'; %factors time series (K x T) --> VAR(1): N=K

T=size(X,2);


mean_X=mean(X,2); %
mean_F=mean(F,2);
 
Xcent=X-repmat(mean_X,1,size(X,2));
Fcent=F-repmat(mean_F,1,size(F,2));
 
cov_XF=(1/T)*(Xcent)*(Fcent)';
cov_F=(1/T)*(Fcent)*(Fcent)'; 

b=cov_XF/cov_F; % loadings
a=mean_X-b*mean_F; %shift (set such that residuals have zero mean)

% residuals
eps=X-repmat(a,1,T)-b*F; %residuals (N x T) 

pick=Mon*(pick_t2m-1)+pick_mon;
IIDAnalysis(eps(pick,:)')

