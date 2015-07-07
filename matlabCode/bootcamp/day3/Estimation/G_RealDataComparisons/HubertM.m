function [Mu,Sigma]=HubertM(x)
% this function computes the location and scatter parameters
% from T N-dimensional joint i.i.d realizations of a multivariate random
% variable as represented by the TxN panel x
% see Sec. 4.5.3 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci


Tolerance=10^(-6);
Error=10^6;
[T,N]=size(x);
w=ones(T,1);
Zeros=zeros(N,1);
Mu=Zeros;
Sigma=Zeros*Zeros';

d_0=sqrt(N)+sqrt(2); % cutoff for Huber-Hampel outlier penalizer

while Error>Tolerance

    Mu_Old=Mu;
    Sigma_Old=Sigma;

    % E step
    Mu=Zeros;
    for t=1:T
        Mu=Mu+w(t)*x(t,:)';
    end
    Mu=Mu/sum(w);

    Sigma=Zeros*Zeros';
    for t=1:T
        Sigma=Sigma+w(t)*w(t)*(x(t,:)'-Mu)*(x(t,:)'-Mu)';
    end
    Sigma=Sigma/(w'*w);

    % M step
    InvS=inv(Sigma);
    d=[];
    for t=1:T
        d=[d; sqrt((x(t,:)'-Mu)'*InvS*(x(t,:)'-Mu))];
    end
    w=OutlierCutoff(d,d_0);

    % convergence
    Error = trace((Sigma-Sigma_Old).^2) + (Mu-Mu_Old)'*(Mu-Mu_Old);

end;