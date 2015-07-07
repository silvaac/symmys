function [Mu,Sigma] = MleRecursionForTFast(x,Nu,Tolerance)
% this function computes recursively the ML estimators of location and scatter 
% of a multivariate Student t distribution with given degrees of freedom
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

[T,N]=size(x);
Ones_N=ones(1,N); % fixed for fast matrix operation
Ones_T=ones(T,1); % fixed for fast matrix operation 

% initialize variables
w=ones(T,1);
Mu=0*Ones_N';
Sigma=0*Ones_N'*Ones_N;

Error=10^6;
% start main loop
while Error>Tolerance
    
    % update
    Mu_Old=Mu;
    Sigma_Old=Sigma;
    
    % Step 1
    W=w*Ones_N;
    Mu=sum(W.*x,1)'/sum(w);
    
    x_c=x-Ones_T*Mu';
    Sigma=(W.*x_c)'*x_c/T;

    % Step 2
    InvS=inv(Sigma);
    Ma2=sum((x_c*InvS).*x_c,2);
    w=(Nu+N)./(Nu+Ma2);

    % convergence
    Error = trace((Sigma-Sigma_Old).^2)/N + (Mu-Mu_Old)'*(Mu-Mu_Old)/N;
    
end