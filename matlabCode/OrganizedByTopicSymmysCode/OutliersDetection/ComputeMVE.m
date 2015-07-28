function [MVE_Location,MVE_Dispersion]=ComputeMVE(Data)
% this function computes the minimum volume ellipsoid for a given time series. 
% The location and scatter parameters that define the ellipsoid are 
% multivariate high-breakdown estimators of location and scatter. 
% See Sec. 4.6.1 of "Risk and Asset Allocation" - Springer (2005), by A. Meucci
% for the theory and the routine implemented below


NumObservations=size(Data,1);
Ones=ones(NumObservations,1);
m=mean(Data)';
S=cov(Data);
det_S_New=0;
w=1/NumObservations*Ones;
KeepLoop=1;
while KeepLoop
  Mahalanobis=[];
  for t=1:NumObservations
    x_t=Data(t,:)';
    Mahalanobis = [Mahalanobis (x_t-m)'*inv(S)*(x_t-m)];
  end
  Update=find(Mahalanobis>1);
  w(Update)=w(Update).*Mahalanobis(Update)';
  m=Data'*w/sum(w);
  S=(Data-Ones*m')'*diag(w)*(Data-Ones*m');
  
  det_S_Old=det_S_New;
  det_S_New=det(S);
  KeepLoop=(det_S_Old/det_S_New < .99999);
end
MVE_Location=m;
MVE_Dispersion=S;

