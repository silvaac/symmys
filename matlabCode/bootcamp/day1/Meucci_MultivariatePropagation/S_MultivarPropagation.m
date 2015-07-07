% ========================= 
%  Multivariate diffusion of location-dispersion ellipsoid, as described in
%  Meucci, A., (2010) "Visualizing Risk: Square-Root Rule, Covariances
%  and Ellipsoids", GARP Risk Professional, February
% 
%  Most recent version of article and code available at http://symmys.com/node/137
%  =========================

clear; clc; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs

% data
load DB_swaps 

% aggregation steps in days
Steps=[1 5 22]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(X(:,1),X(:,2),'.');
T=size(X,1);
for s=1:length(Steps)
    
    % compute series at aggregated time steps
    k=Steps(s);
    AggX=[];
    t=1;
    while t+k+1<=T
        NewTerm=sum(X([t : t+k-1],:),1);
        AggX=[AggX
            NewTerm];
        t=t+k;
    end
    
    % empirical mean/covariance
    Agg(s).M_hat=mean(AggX)';
    Agg(s).S_hat=cov(AggX);
    
    % mean/covariance implied by propagation law of risk for invariants
    Agg(s).M_norm=k/Steps(1)*Agg(1).M_hat;
    Agg(s).S_norm=k/Steps(1)*Agg(1).S_hat;

    % plots
    hold on
    h1=TwoDimEllipsoid(Agg(s).M_norm,Agg(s).S_norm,1,0,0);
    set(h1,'color','k','linewidth',1,'linestyle','--')
    
    hold on
    h2=TwoDimEllipsoid(Agg(s).M_hat,Agg(s).S_hat,1,0,0);
    set(h2,'color','r','linewidth',2)
    
end
xlabel( Names{1} )
ylabel( Names{2} )
grid off
