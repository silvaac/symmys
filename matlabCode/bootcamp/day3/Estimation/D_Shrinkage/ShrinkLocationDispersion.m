function [Mu_shr, Sigma_shr]=ShrinkLocationDispersion(X)

[T,N] = size(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% standard sample estimate
Mu_hat=mean(X)';
X = X - ones(T,1) * Mu_hat';
Sigma_hat = X'*X/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shrinkage of location parameter to fixed target 
m=0*ones(N,1); % fixed target 
% compute optimal weight
Lambda_hat=eig(Sigma_hat);
a=1/T*(sum(Lambda_hat)-2*max(Lambda_hat))/((Mu_hat-m)'*(Mu_hat-m)); 
a=max(0,min(a,1));     % restrict to sensible weights

Mu_shr=(1-a)*Mu_hat+a*m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shrinkage of dispersion parameter to sphere
C=mean(Lambda_hat)*eye(N);
Numerator=0;
for t=1:T
    Numerator = Numerator +  1/T*trace(  (X(t,:)'*X(t,:)-Sigma_hat)^2  ) ;
end
Denominator=trace( (Sigma_hat-C)^2);
a=1/T*Numerator/Denominator;
a=max(0,min(a,1)); % restrict to sensible weights

Sigma_shr=(1-a)*Sigma_hat+a*C;