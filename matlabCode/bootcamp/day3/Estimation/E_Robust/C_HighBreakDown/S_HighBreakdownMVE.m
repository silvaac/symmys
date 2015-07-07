% this script shows how to use the minimum volume ellipsoid to detect outliers. 
% See Sec. 4.6.1 of "Risk and Asset Allocation" - Springer (2005), by A. Meucci
% for the theory and the routine implemented below

clear;  close all;  clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate observations
Mu=[0; 0];
r=-.90;
sig=[1 1];
T=50;

Outliers=10*rand(7,2);

Sigma=diag(sig)*[1 r; r 1]*diag(sig);
% "good" observations
Sample = mvnrnd(Mu,Sigma,T);
% add "bad" observation(s)
CorruptSample=[Sample
    Outliers];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute high-breakdown estimates
T=size(CorruptSample,1);
Index=[1:T];

Sample_Length=[];
Vol_MVE=[];
for j=1:ceil(T/2)
    
    % compute MVE
    [MVE_Location,MVE_Dispersion]=ComputeMVE(CorruptSample);

    % store results
    Store(j).MVE_Location=MVE_Location;
    Store(j).MVE_Dispersion=MVE_Dispersion;
    
    Sample_Length=[Sample_Length size(CorruptSample,1)];
    Vol_MVE=[Vol_MVE sqrt(det(MVE_Dispersion))];
    
    Store(j).Index=Index;

   % erase one outlier
    Rejected=RejectOutlier(CorruptSample,Index);
    CorruptSample(Rejected,:)=[];
    Index(Rejected)=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results

figure % data
h=plot(Sample(:,1),Sample(:,2),'.');
set(h,'color','k')
hold on
h=plot(Outliers(:,1),Outliers(:,2),'.');
set(h,'color','r')
for j=1:ceil(T/2)
    hold on
    E=TwoDimEllipsoid(Store(j).MVE_Location,Store(j).MVE_Dispersion,1,0,0);
    set(E,'color','b','linewidth',1)
end
title('MVE estimates')

figure % volume of ellipsoid as function of sample length
plot(Sample_Length,Vol_MVE)
grid on 