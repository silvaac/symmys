% this script computes the PCA structure of the swap markets in the continuum limit
% see Sec. 3.5. in "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MinMaturity=3;   % select section of the zero-swap curve for the empirical analysis (min 0 max 10)
MaxMaturity=10;

Frequency_Range=[0 : 1/252 : .2]; % select frequency range for the display of the contiuous limit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load DB_SwapCurve
SelectCurve=find(Maturities<=MaxMaturity & Maturities>=MinMaturity );
Zero_Spot=-log(DF(:,SelectCurve))./(ones(length(Dates),1)*Maturities(SelectCurve));
Changes_Zero_Spot=Zero_Spot(2:end,:)-Zero_Spot(1:end-1,:);
Covariance=Changes_Zero_Spot'*Changes_Zero_Spot/size(Changes_Zero_Spot,1);
[StDeviations,Correlation]=cov2corr(Covariance);

[X,Y]=meshgrid(Maturities(SelectCurve),Maturities(SelectCurve));
x0=.1;
options=[];
Gamma = fminunc('FitError',x0,options,X,Y,Correlation);
Correlation_fit=exp(-Gamma*abs(X-Y));

% plot emprical and continuum correlation 
X_Lim=[min(min(X)) max(max(X))];
Y_Lim=[min(min(Y)) max(max(Y))];
Z_Lim=[.85 1.05];

figure
subplot('Position',[0.05 .55 .9 .40] )
m=min(min(Correlation)); M=max(max(Correlation));
uu= (Correlation-m)/(M-m);
hl=surf(X,Y,Correlation);
set(hl,'CData',uu,'EdgeColor',[.5 .5 .5],'meshstyle','row')
set(gca,'xlim',X_Lim,'ylim',Y_Lim,'zlim',Z_Lim);
caxis([-.1 .9])
grid on
box on

subplot('Position',[0.05 .05 .9 .40] )
m=min(min(Correlation_fit)); M=max(max(Correlation_fit));
uu= (Correlation_fit-m)/(M-m);
hl=surf(X,Y,Correlation_fit);
set(hl,'CData',uu,'EdgeColor',[.5 .5 .5],'meshstyle','row')
grid off
set(gca,'xlim',X_Lim,'ylim',Y_Lim,'zlim',Z_Lim);
caxis([-.1 .9])
grid on
box on

colormap('gray')

% plot normalized eigenvalues, r-square and eigenfuctions
figure
subplot('Position',[.05 .75 .9 .24])
EigValues=Gamma*Gamma./(Gamma*Gamma+Frequency_Range.*Frequency_Range);
x=[Frequency_Range 0];
y=[EigValues 0];
fill(x,y,[.7 .7 .7])
hold on 
h2=plot(Frequency_Range,EigValues);
set(h2,'linewidth',2,'color','k');
set(gca,'xlim',[Frequency_Range(1) Frequency_Range(end)],'ylim',[-.1 1.1]);
grid off
box off

subplot('Position',[.05 .4 .9 .24])
R_2=2/pi*atan(Frequency_Range/Gamma);
x=[Frequency_Range Frequency_Range(end)];
y=[R_2 0];
fill(x,y,[.7 .7 .7])
hold on 
h2=plot(Frequency_Range,R_2);
set(h2,'linewidth',2,'color','k');
set(gca,'xlim',[Frequency_Range(1) Frequency_Range(end)],'ylim',[-.1 1.1]);
grid off
box off


subplot('Position',[.05 .05 .9 .24])
CutOff=.2;
Steps=3;
MaturityDiff=[0:30];
Omegas=[0: CutOff/(Steps-1) : CutOff];
maxY=0; minY=0;
for i=1:length(Omegas)
    a=cos(Omegas(i)*MaturityDiff);
    EigFunction=a/sqrt(a*a');
    maxY=max(maxY,max(EigFunction));
    minY=min(minY,min(EigFunction));
    hold on 
    h=plot(abs(MaturityDiff),EigFunction);
    set(h,'linewidth',2,'color','k')
end
set(gca,'ylim',[minY-.1*(maxY-minY) maxY+.1*(maxY-minY)]);
grid off
box off