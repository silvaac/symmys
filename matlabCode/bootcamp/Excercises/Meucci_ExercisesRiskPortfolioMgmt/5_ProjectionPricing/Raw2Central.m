function mu=Raw2Central(mu_)

% code to map raw moments into central moments 
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci and 
% "Exercises in advanced risk and portfolio management - with step-by-step solutions and fully documented code", 
% free E-Book available at http://ssrn.com/abstract=1447443

% note: first central moment defined as first non-central moment, i.e. the expectation

N=length(mu_);
mu=mu_;

for n=2:N
    mu(n) = ((-1)^n) * (mu_(1))^(n);
    for k=1:n-1
        mu(n) =  mu(n) + nchoosek(n,k) * ((-1)^(n-k)) * mu_(k) * (mu_(1))^(n-k);
    end
    mu(n) = mu(n)+ mu_(n);
end