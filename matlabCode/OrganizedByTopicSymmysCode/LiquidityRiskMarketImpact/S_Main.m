%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script computes the liquidity-risk and funding-risk adjusted P&L distribution 
% see A. Meucci - "A Fully Integrated Liquidity and Market Risk Model", Financial Analyst Journal, 68, 6, 35-47 (2012)
% Code by A. Meucci. This version August 2013. 

% Last version of code and article available at http://symmys.com/node/350
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; clc;  close all;

% INPUTS
%*********************************************************************
% liquidation policy at horizon as fraction of investment
Policy=-1;

% collinearity of liquidity perturbations
CollinLiq=1;

% select only some stock in portfolio and equally allocate capital as fraction of daily dollar volume
Selectstock=[1:10]';
Capital_perDailyVolume= .2;


% PREPARE DATA
%*********************************************************************
% load  Daily_Prices: closing prices 
%       Daily_Volumes_Shares: daily volumes   
%       Daily_Liq: Morgan Stanley liquidity index 
load DB

% Prices and returns
%Daily_Prices = Daily_Prices(:,Selectstock);
Prices_0 = Daily_Prices(end,:)';
Daily_LogRets=log(Daily_Prices(1:end-1,:)./Daily_Prices(2:end,:));
[J N]=size(Daily_LogRets);

% volumes in shares
%Daily_Volumes = Daily_Volumes_Shares(:,Selectstock);
Volumes_0 =  Daily_Volumes_Shares(end,:)';
Volumes_t =mean(Daily_Volumes_Shares(end-250:end,:))';

% liquidity index
Daily_LiqChanges=diff(Daily_Liq);
Liq_0 = Daily_Liq(end,:)';

% normal simulations
    X=[Daily_LogRets  Daily_LiqChanges];
    m_X=mean(X);
    s2_X=cov(X); %covariance
    J=100000;
    RandStream.setGlobalStream(RandStream('mt19937ar', 'seed', 11)); 
    X=mvnrnd(m_X,s2_X,J);  
    Daily_LogRets=X(:,1:N);
    Daily_LiqChanges=X(:,end);

% Fully Flexible Probabilties associated with each scenario
Probs=ones(J,1)/J;

% stock prices at horizon 
Prices_t = repmat(Prices_0',J,1).*exp(Daily_LogRets);

% liquidity index at horizon
Liq_t = Liq_0 .*exp(Daily_LiqChanges);

% pure market risk: p&L due to market risk
PnL_mkt = Prices_t - repmat(Prices_0',J,1);


% PORTFOLIO COMPUTATIONS
%*********************************************************************
% portfolio and liquidation policy
Weights=zeros(N,1);
Weights(Selectstock)=1/length(Selectstock);
DollarVolume_0 = Volumes_0'*Prices_0;
Capital=Capital_perDailyVolume*DollarVolume_0;

h=Capital*Weights./Prices_0;

PnL_mkt_h=PnL_mkt*h;                           

% LIQUIDITY ADJUSTMENT 
%*********************************************************************
% liquidation policy
Dh=Policy*h;

% market impact
b_a=0.01*ones(N,1);                          
Linear=-b_a.*Prices_0.*abs(Dh);
NonLinear=-(10^5)*Prices_0.*std(Daily_LogRets)'.*((abs(Dh)./Volumes_t ).^1.5);         
m_Dh = Linear + NonLinear;                                      
    
% state-dependent expected liquidity impact on all stocks
s_g1=.5*std(PnL_mkt_h);
g1=-min(PnL_mkt_h,-s_g1)/s_g1;
m_Dh_x=repmat(g1,1,N).*repmat(m_Dh',J,1);    % (14)

% state-dependent expected liquidity impact on portfolio
m_Dh_h=m_Dh_x*ones(N,1);                    % (23)

% state-independent uncertainty on liquidity impact on portfolio
s_Dh=1.5*m_Dh;                                                % 
r2_Dh=(1-CollinLiq)*corr(Daily_LogRets)+CollinLiq*ones(N,N);  % 
s2_Dh=diag(s_Dh)*r2_Dh*diag(s_Dh);                            % 
s2_Dh_h=ones(N,1)'*s2_Dh*ones(N,1);                           % 
s_Dh_h=sqrt(s2_Dh_h);
s_Dh_h=max(s_Dh_h,.01*std(PnL_mkt_h));                        % regularization

% TOTAL P&L
%*********************************************************************
% conditional center and scatter
m_j=PnL_mkt_h+m_Dh_h;
s_j=s_Dh_h*ones(J,1);

% pdf and cdf: taking and not taking into account funding cost
nu=100;
f_Pi = @(x) (Probs./s_j)'*tpdf((x-m_j)./s_j,nu);
F_Pi = @(x) Probs'*tcdf((x-m_j)./s_j,nu);

NGrid=200;
x_=linspace(min(PnL_mkt_h)-s_Dh_h,max(PnL_mkt_h)+s_Dh_h, NGrid);
p_=[];
f_Pi_plot=[];
f_Pi_funding_plot=[];
for k=1:NGrid
    p_=[p_
        F_Pi(x_(k))];
    f_Pi_plot=[f_Pi_plot
        f_Pi(x_(k))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure
NumBins=round(10*log(J));
[n,Ds]=hist(PnL_mkt_h,NumBins); % compute bin width
D=Ds(2)-Ds(1);
hh=bar(Ds,n/(J*D)); % plot histogram
set(hh,'edgecolor',.4*[1 1 1],'facecolor',.4*[1 1 1])
hold on
hh1=plot(x_,f_Pi_plot,'r');
legend('pure market P&L','market + liquidity P&L')
xlim([-10^9 10^9])