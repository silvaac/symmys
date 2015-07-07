function [ga,mu]=SummStats(X,N)

% compute central moments
mu=zeros(1,N);
mu(1)=mean(X);
for n=2:N
    mu(n)=moment(X,n);
end

% compute standardized statistics 
ga=mu;
ga(2)=sqrt(mu(2));
for n=3:N
    ga(n)=mu(n)/(ga(2)^n);
end


