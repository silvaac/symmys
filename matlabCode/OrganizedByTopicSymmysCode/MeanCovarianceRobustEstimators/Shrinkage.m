function [M_Shr,S_Shr]=Shrinkage(X)

[T,N]=size(X);

M_NoPar=mean(X)';
S_NoPar=cov(X);

b=0*ones(N,1); % target 
Lambda_hat=eig(S_NoPar); % compute optimal weight
a=1/T*(sum(Lambda_hat)-2*max(Lambda_hat))/((M_NoPar-b)'*(M_NoPar-b)); 
a=max(0,min(a,1));   % restrict to sensible weight
M_Shr=(1-a)*M_NoPar+a*b;

C=mean(Lambda_hat)*eye(N); % target
% compute optimal weight
Numerator=0;
for t=1:T
    Numerator = Numerator +  1/T*trace(  (X(t,:)'*X(t,:)-S_NoPar)^2  ) ;
end
Denominator=trace( (S_NoPar-C)^2);
a=1/T*Numerator/Denominator;
a=max(0,min(a,1)); % restrict to sensible weight
S_Shr=(1-a)*S_NoPar+a*C;