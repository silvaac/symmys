% this script performs the principal component analysis of the swap curve. it computes and plots, among others, 
% 1. the location-dispersion ellipsoid of three rates along with the 3-d location-dispersion ellipsoid
% 2. the spectrum profile, the profile of the r-square, and the Fourier-basis profile of the eigenvectors
% 3. the empirical distribtuion of the first three uncorrelated principal factors
% 4. the effect on the curve of the first three uncorrelated principal factors 
% see Sec. 3.5. in "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clear; clc; close all;
load DB_SwapCurve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs 
MinMaturity=2;   % min = 0 max=10
MaxMaturity=10;
Nodes=1;         % (.25=1 quarter,  1=1 year,  2=2 years)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computations
SelectCurve=find(Maturities<=MaxMaturity & Maturities>=MinMaturity );
SelectCurve = intersect(SelectCurve,[0:4*Nodes:40]);
Maturities=Maturities(SelectCurve);

Zero_Spot=-log(DF(:,SelectCurve))./(ones(length(Dates),1)*Maturities);
Current_Zero_Spot=Zero_Spot(end,:);

Changes_Zero_Spot=Zero_Spot(1:end-1,:)-Zero_Spot(2:end,:);
Covariance=Changes_Zero_Spot'*Changes_Zero_Spot/size(Changes_Zero_Spot,1); % estimate covariance by assuming mean is zero
[StDeviations,Correlation]=cov2corr(Covariance);

% perform principal component analysis
[EigVectors,EigValues]=eig(Covariance);
[dummy,Index]=sort(-diag(EigValues));
EigVectors=EigVectors(:,Index);
EigValues=diag(EigValues(Index,Index));

EigV1=EigVectors(:,1);
EigV2=EigVectors(:,2);
EigV3=EigVectors(:,3);
Factor1=Changes_Zero_Spot*EigV1;
Factor2=Changes_Zero_Spot*EigV2;
Factor3=Changes_Zero_Spot*EigV3;
Std1=sqrt(EigValues(1));
Std2=sqrt(EigValues(2));
Std3=sqrt(EigValues(3));

R_Square=[];
for i=1:length(EigValues);
    R_Square=[R_Square sum(EigValues(1:i))];
end
R_Square=R_Square/sum(EigValues);

% select data for 3-d scatter plot and compute location-dispersion ellipsoid
L=size(Changes_Zero_Spot,2);
Y=10000*Changes_Zero_Spot(:,[1 round(L/2) L]);
S=cov(Y);
[EVrs,EVls] = pcacov(S);
Theta = [-pi/2: pi/50 : 0];
Phi = [0 : pi/80 : 2*pi];
Scale=3;
Ellipsoid_x=zeros(length(Theta),length(Phi)); Ellipsoid_y=Ellipsoid_x; Ellipsoid_z=Ellipsoid_x;
for i=1:length(Theta)
    for j=1:length(Phi)
        
        y=[cos(Theta(i));  sin(Theta(i))*sin(Phi(j)); sin(Theta(i))*cos(Phi(j))];
        Ellipsoid = Scale*EVrs*diag(sqrt(EVls))*y;
        
        Ellipsoid_x(i,j)= Ellipsoid(1);
        Ellipsoid_y(i,j)= Ellipsoid(2);
        Ellipsoid_z(i,j)= Ellipsoid(3);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

% plot covariance matrix
figure
[m1,m2]=meshgrid(Maturities,Maturities);
[aa,C]=cov2corr(Covariance);
surf(m1,m2,Covariance)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot rate changes and location-dispersion ellipsoid
figure
h=plot3(Y(:,1),Y(:,2),Y(:,3),'.')
set(h,'color','k')
hold on
m=min(min(Ellipsoid_z)); M=max(max(Ellipsoid_z));
uu= (Ellipsoid_z-m)/(M-m);
hh=surf(Ellipsoid_x,Ellipsoid_y,Ellipsoid_z);
set(hh,'CData',uu,'EdgeColor',[.6 .6 .6],'meshstyle','row')
set(gca,'xdir','reverse','ydir','reverse')
ax1=[0 0 0; 0 35 0];
ax2=[0 0 0; 0 0 35];
ax3=[0 0 0; 35 0 0];
h=plot3(ax1(:,1),ax1(:,2),ax1(:,3));
set(h,'color',[.3 .3 .3],'linewidth',2)
h=plot3(ax2(:,1),ax2(:,2),ax2(:,3));
set(h,'color',[.3 .3 .3],'linewidth',2)
h=plot3(ax3(:,1),ax3(:,2),ax3(:,3));
set(h,'color',[.3 .3 .3],'linewidth',2)
axis equal
set(gca,'xlim',[-40 40],'ylim',[-40 40],'zlim',[-40 40],'Position',[-.1 -.1  1.2  1.2]);
grid on
box on
axis xy
colormap('gray')
caxis([-1.5 1.1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot normalized eigenvalues, r-square and eigenfuctions
figure 
subplot('Position',[.05 .75 .9 .24])
h=bar(EigValues/EigValues(1));
set(h,'FaceColor',[.7 .7 .7],'EdgeColor','k')
set(gca,'ytick',[],'ylim',[0 1.1]);
grid off
box off

subplot('Position',[.05 .4 .9 .24])
h=bar(R_Square);
set(h,'FaceColor',[.7 .7 .7],'EdgeColor','k')
set(gca,'ylim',[.9 1.01]);
grid off
box off

subplot('Position',[.05 .05 .9 .24])
h=plot(Maturities,[EigV1 EigV2 EigV3]);
set(h,'linewidth',2,'color','k')
minY=min(min([EigV1 EigV2 EigV3])); maxY=max(max([EigV1 EigV2 EigV3]));
set(gca,'ylim',[minY-.1*(maxY-minY) maxY+.1*(maxY-minY)]);
grid off
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot factors empirical distribution
figure
Nbins=round(9*log(length(Factor1)));

subplot('Position',[.05 .75 .9 .24])
h=histfit(10000*Factor1,Nbins);
Xs=get(h(1),'XData'); Ys=get(h(1),'YData');
Ys=Ys/(max(max(abs(Ys))));
set(h(1),'XData',Xs,'YData',Ys,'FaceColor',[.9 .9 .9],'EdgeColor','k')
Xs=get(h(2),'XData'); Ys=get(h(2),'YData');
Ys=Ys/(max(max(abs(Ys))));
set(h(2),'XData',Xs,'YData',Ys,'Color','k','linewidth',2)
set(gca,'xlim',[-150 150],'ylim',[-.01 1.1]);
grid off
box off

subplot('Position',[.05 .41 .9 .24])
h=histfit(10000*Factor2,Nbins);
Xs=get(h(1),'XData'); Ys=get(h(1),'YData');
Ys=Ys/(max(max(abs(Ys))));
set(h(1),'XData',Xs,'YData',Ys,'FaceColor',[.9 .9 .9],'EdgeColor','k')
Xs=get(h(2),'XData'); Ys=get(h(2),'YData');
Ys=Ys/(max(max(abs(Ys))));
set(h(2),'XData',Xs,'YData',Ys,'Color','k','linewidth',2)
set(gca,'xlim',[-40 40],'ylim',[-.01 1.1]);
grid off
box off

subplot('Position',[.05 .07 .9 .24])
h=histfit(10000*Factor3,Nbins);
Xs=get(h(1),'XData'); Ys=get(h(1),'YData');
Ys=Ys/(max(max(abs(Ys))));
set(h(1),'XData',Xs,'YData',Ys,'FaceColor',[.9 .9 .9],'EdgeColor','k')
Xs=get(h(2),'XData'); Ys=get(h(2),'YData');
Ys=Ys/(max(max(abs(Ys))));
set(h(2),'XData',Xs,'YData',Ys,'Color','k','linewidth',2)
set(gca,'xlim',[-15 15],'ylim',[-.01 1.1]);
grid off
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot 3-std effect of each factor on the curve 
figure 
subplot('Position',[.05 .75 .9 .24])
Shift=3*Std1*EigV1';
h=plot(Maturities,100*Current_Zero_Spot);
set(h,'linewidth',1.5,'color','k')
hold on 
h=plot(Maturities,100*(Current_Zero_Spot+Shift));
set(h,'linewidth',1.5,'color','k')
hold on 
h=plot(Maturities,100*(Current_Zero_Spot-Shift));
set(h,'linewidth',1.5,'color','k')
set(gca,'ylim',[7.5 9.5]);
grid on


subplot('Position',[.05 .41 .9 .24])
Slope=3*Std2*EigV2';
h=plot(Maturities,100*Current_Zero_Spot);
set(h,'linewidth',1.5,'color','k')
hold on 
h=plot(Maturities,100*(Current_Zero_Spot+Slope));
set(h,'linewidth',1.5,'color','k')
hold on 
h=plot(Maturities,100*(Current_Zero_Spot-Slope));
set(h,'linewidth',1.5,'color','k')
set(gca,'ylim',[8 9]);
grid on


subplot('Position',[.05 .07 .9 .24])
Hump=3*Std3*EigV3';
h=plot(Maturities,100*Current_Zero_Spot);
set(h,'linewidth',1.5,'color','k')
hold on 
h=plot(Maturities,100*(Current_Zero_Spot+Hump));
set(h,'linewidth',1.5,'color','k')
hold on 
h=plot(Maturities,100*(Current_Zero_Spot-Hump));
set(h,'linewidth',1.5,'color','k')
set(gca,'ylim',[8.2 8.85]);
grid on