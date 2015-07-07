function Rejected=RejectOutlier(Sample,Index);
% this function finds the "worst" outlier in a time series. 
% See Sec. 4.6.1 of "Risk and Asset Allocation" - Springer (2005), by A. Meucci
% for the theory and the routine implemented below

T=size(Sample,1);
m=mean(Sample)';
U=Sample-ones(T,1)*m';
Lambda=diag(U*inv(U'*U)*U');
[a,Rejected]=max(Lambda);
