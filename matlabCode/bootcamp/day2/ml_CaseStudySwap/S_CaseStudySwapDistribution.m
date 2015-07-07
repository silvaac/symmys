% this script computes the PnL distribution at the investment horizon of a forward swap position
% three methods are compared: 
% 1. Monte Carlo exact pricing of the distribution of the first three PCA factors
% 2. First-order (dv01) approximation 
% 3. Second-order (convexity) approximation 
% see Sec. 3.5. in "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clear; clc; close all;
load DB_SwapCurve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs 
MinMaturity=2;  MaxMaturity=10;  % min = 1 max=10

Fixed_Rate=.049;
Horizon=1/12;
Notional=1000000;
NumSimul=100000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computations

% run statistical analysis on year-apart rates
Estimation_Interval=1/52;  % cannot change this!!!
Time2Mats=[MinMaturity : .25 : MaxMaturity];
Time2Mats_Hor=Time2Mats-Horizon;

SelectCurve=find(Maturities<=MaxMaturity & Maturities>=MinMaturity);  
Maturities=Maturities(SelectCurve);
Zero_Spot=-log(DF(:,SelectCurve))./(ones(length(Dates),1)*Maturities);
InterpZero_Spot=interp1(Maturities,Zero_Spot',Time2Mats_Hor,'spline','extrap')';
Changes_Zero_Spot=InterpZero_Spot(1:end-1,:)-InterpZero_Spot(2:end,:);
Covariance=cov(Changes_Zero_Spot);
Mean=mean(Changes_Zero_Spot);
[StDeviations,Correlation]=cov2corr(Covariance);
[EigVectors,EigValues]=pcacov(Covariance);

EigV1=-EigVectors(:,1); 
Ones=ones(length(EigV1),1);
Alpha=EigV1'*Ones/(EigV1'*EigV1);
Var_1=EigValues(1)/(Alpha*Alpha);
Factor_Loadings=[Ones EigVectors(:,2) EigVectors(:,3)];
Factor_Covariance=diag([Var_1  EigValues(2) EigValues(3)]);

Coefficients=Fixed_Rate*.25*Ones;
Coefficients(end)=Coefficients(end)+1;
Coefficients(1)=-1;

Current_ZeroSpot=interp1(Maturities,Zero_Spot(1,:),Time2Mats,'spline','extrap');
Bond_Prices=exp(-Current_ZeroSpot.*Time2Mats)*Notional;
Current_ZeroSpot_Hor=InterpZero_Spot(1,:);
Bond_Prices_Slide=exp(-Current_ZeroSpot_Hor.*Time2Mats_Hor)*Notional;

Swap_CurrentValue=Bond_Prices*Coefficients;
Slide=Bond_Prices_Slide*Coefficients;
PVBP=(Bond_Prices_Slide.*Time2Mats_Hor)*Coefficients;
Convexity=.5*(Bond_Prices_Slide.*Time2Mats_Hor.*Time2Mats_Hor)*Coefficients;

% 3 factors simulation
Mu=0*Ones;
Sigma=Factor_Loadings*Factor_Covariance*Factor_Loadings'*Horizon/Estimation_Interval;
X= mvnrnd(Mu,Sigma,NumSimul/2);
X=[X
   -X];
Swap_HorizonValue=exp(X.*(ones(NumSimul,1)*Time2Mats_Hor))*(Coefficients'.*Bond_Prices)';

% order 1 approximation
BottomRange_P=-Notional/30*sqrt(Horizon/Estimation_Interval);
TopRange_P=Notional/30*sqrt(Horizon/Estimation_Interval);
Step=(TopRange_P-BottomRange_P)/1000;
PriceRange=[BottomRange_P : Step : TopRange_P];
PriceDens_OrdOne=normpdf(PriceRange,Slide,sqrt(PVBP^2*Var_1*Horizon/Estimation_Interval));

% order 2 approximation
g=(Slide-PVBP*PVBP/(4*Convexity));
q=Convexity*Var_1*Horizon/Estimation_Interval;
DegFreedom=1;
Non_Centrality=(  -PVBP/(2*Convexity*sqrt(Var_1*Horizon/Estimation_Interval))  )^2;

PriceDens_OrdTwo=1/q*ncx2pdf((PriceRange-g)/q,DegFreedom,Non_Centrality);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

figure

subplot('Position',[.05 .75 .9 .24]) % 1-factor, first-order approximation
y=PriceDens_OrdOne/max(PriceDens_OrdOne);
fill(PriceRange,y,[.7 .7 .7])
hold on 
h=plot(PriceRange,y);
set(h,'color','k','linewidth',2)
hold on
h=plot(Swap_CurrentValue,0,'d');   % plot current value
set(h,'markersize',5,'markerfacecolor','k','markeredgecolor','k')
hold on
h=plot(Slide,0,'s');               % plot slide
set(h,'markersize',5,'markerfacecolor','k','markeredgecolor','k')
set(gca,'ytick',[],'xlim',[PriceRange(1) PriceRange(end)],'ylim',[-.1 1.1]);
grid on

subplot('Position',[.05 .41 .9 .24]) % 1-factor, second-order approximation
y=PriceDens_OrdTwo/max(PriceDens_OrdTwo);
fill(PriceRange,y,[.7 .7 .7])
hold on
h=plot(PriceRange,y);
set(h,'color','k','linewidth',2)
hold on
h=plot(Swap_CurrentValue,0,'d');   % plot current value
set(h,'markersize',5,'markerfacecolor','k','markeredgecolor','k')
hold on
h=plot(Slide,0,'s');               % plot slide
set(h,'markersize',5,'markerfacecolor','k','markeredgecolor','k')
set(gca,'ytick',[],'xlim',[PriceRange(1) PriceRange(end)],'ylim',[-.1 1.1]);
grid on

subplot('Position',[.05 .07 .9 .24]) % 3-factor simulation
NumBins=round(10*log(NumSimul));
[n,xout] = hist(Swap_HorizonValue,NumBins) ;
h=bar(xout,n,1);
Xs=get(h,'XData'); Ys=get(h,'YData');
Ys=Ys/(max(max(abs(Ys))));
set(h,'XData',Xs,'YData',Ys,'FaceColor',[.7 .7 .7],'EdgeColor','k')
hold on
h=plot(Swap_CurrentValue,0,'d');   % plot current value
set(h,'markersize',5,'markerfacecolor','k','markeredgecolor','k')
hold on
h=plot(Slide,0,'s');               % plot slide
set(h,'markersize',5,'markerfacecolor','k','markeredgecolor','k')
set(gca,'ytick',[],'xlim',[PriceRange(1) PriceRange(end)],'ylim',[-.1 1.1]);
grid on