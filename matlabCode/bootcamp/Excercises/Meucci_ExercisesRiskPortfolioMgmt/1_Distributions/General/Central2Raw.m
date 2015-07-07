function mu_=Central2Raw(mu)

% code to map central moments into raw moments
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci and 
% "Exercises in advanced risk and portfolio management - with step-by-step solutions and fully documented code", 
% free E-Book available at http://ssrn.com/abstract=1447443

% note: first central moment defined as first non-central moment, i.e. the expectation

N=length(mu);
mu_=mu;

for n=2:N
    mu_(n) = ((-1)^(n+1)) * (mu(1))^(n);
    for k=1:n-1
        mu_(n) =  mu_(n) + nchoosek(n,k) * ((-1)^(n-k+1)) * mu_(k) * (mu_(1))^(n-k);
    end
    mu_(n) = mu_(n)+ mu(n);
end