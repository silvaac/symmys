clc; clear; close all;

load db_BondAttribution 
% B = key rate durations
% F = key rate weekly changes
% X = bonds returns net of carry

[T,K,N]=size(B);

U=0*X;
for t=1:T
    U(t,:)=X(t,:)-F(t,:)*squeeze(B(t,:,:));
end

C=corr([U F]);

C_U=C(1:N,1:N)
C_FU=C(1:N,N+1:end)