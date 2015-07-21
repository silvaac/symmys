function Mom = MVOU_Prior(t, x0, theta, sig2, mu)

% This function computes the conditional mean and covariance matrix at different
% monitoring dates t1,t2,...,t_ of  the process
% X_{t1,t2...,t_} =(X_t1, 
%                   X_t2,
%                   .
%                   .
%                   X_t_)
% 
% X_t follows a MVOU process: dX_t = (-theta*X_t+mu)dt + sig*dB_t

%INPUT
% t = vector of monitoring dates.                   [t_ x 1]
% x0 = observation at time 0.                       [n_ x 1]
% theta = transition matrix.                        [n_ x n_]
% sig2 = covariance matrix.                         [n_ x n_]
% mu = vector of drift parameters                   [n_ x 1]

%OUTPUT
%Mom.monitoring_time = monitoring times             [t_*n_ x 1]
%Mom.dimension = labels of the risk drivers         [t_*n_ x 1]    
%Mom.cov = covariance matrix of X_{t1,t2...t_}      [t_*n_ x t_*n_]
%Mom.mean = vector of the means of X_{t1,t2...t_}   [t_*n_ x 1]
%Mom.mean_cost                                      [t_*n_ x 1]
%Mom.mean_lin                                       [t_*n_ x 1]
%Mom.mean_cost and Mom.mean_lin are such that Mom.mean_cost + Mom.mean_lin*x0  = Mom.mean        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tol_eigb = 10^-8; 
t_ = length(t);
n_ = length(x0);
[T,N] = meshgrid(t,1:n_);
Mom = struct('monitoring_time', T(:), 'dimension',N(:),'mean', NaN(t_*n_,1), 'cov',NaN(n_*t_,n_*t_),'mean_cost',NaN(t_*n_,1),'mean_lin',NaN(t_*n_,n_));
kronsum = kron(theta,eye(n_)) + kron(eye(n_),theta);    
[V, D] = eig(kronsum);
lambda = diag(D);
lambda_A = NaN(length(D),1);

M = NaN(n_,t_);
mean_lin = NaN(n_,n_,t_);
mean_cost = NaN(n_,t_);
[V1,D1] = eig(theta);
theta_diag = diag(D1);
F = NaN(n_,1);

for i = 1:t_
    F(theta_diag<=Tol_eigb) = t(i);
    F(theta_diag>Tol_eigb) = (1-exp(-theta_diag(theta_diag>Tol_eigb)*t(i)))./theta_diag(theta_diag>Tol_eigb);
    E = expm(-theta*t(i));
    M(1:n_,i) = (E*x0+ V1*diag(F)*pinv(V1)*mu);        
    M(:,i) = real(M(:,i));
    
    mean_lin(1:n_,1:n_,i) = E';        
    mean_lin(:,:,i) = real(mean_lin(:,:,i));
    mean_cost(1:n_,i) = (V1*diag(F)*pinv(V1)*mu);        
    mean_cost(:,i) = real(mean_cost(:,i));

    vecsig2 = reshape(sig2, n_^2,1);
    lambda_A((abs(lambda) <= Tol_eigb)) = t(i);
    index = abs(lambda) > Tol_eigb;
    lambda_A(index) = (1-exp(-lambda(index)*t(i)))./lambda(index);    
    A = V*diag(lambda_A)/V;
    vecsig2_t = A*vecsig2;
    sig2_t = reshape(vecsig2_t,n_,n_);                              
    sig2_t = real(sig2_t);

    for j = i:t_    
        Mom.cov((i-1)*n_+(1:n_),(j-1)*n_+(1:n_)) = sig2_t*expm(-theta'*(t(j)-t(i)));        
        Mom.cov((j-1)*n_+(1:n_),(i-1)*n_+(1:n_)) = expm(-theta*(t(j)-t(i)))*sig2_t;        
    end    
end

%Note that Mom.mean = Mom.mean_cost + Mom.mean_lin*x0 
Mom.mean = M(:);
Mom.mean_lin = mean_lin(:,:)';
Mom.mean_cost = mean_cost(:);




