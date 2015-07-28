function Posterior = MVOU_Posterior(Mom, Views)

%This function computes the posterior conditional expectation and covariance matrix of
% X_{t1,t2...,t_} =(X_t1, 
%                   X_t2,
%                   .
%                   .
%                   X_t_)
% where t1,t2...,t_ are monitoring times and X_{t1,t2...,t_} is multivariate
% normal distributed.
% The Views are: E{v_mu*X}=mu_view and Cov{v_sig*X}=sig2_view.  

%INPUT
% Mom = struct('monitoring_time','dimension','mean', 'cov','mean_cost','mean_lin');
% where:
% Mom.monitoring_time = monitoring times                        [t_*n_ x 1]
% Mom.dimension = labels of the risk drivers                    [t_*n_ x 1]    
% Mom.cov = prior covariance matrix of X_{t1,t2...t_}           [t_*n_ x t_*n_]
% Mom.mean = prior vector of the means of X_{t1,t2...t_}        [t_*n_ x 1]
% Mom.mean_cost                                                 [t_*n_ x 1]
% Mom.mean_lin                                                  [t_*n_ x 1]
% Mom.mean_cost and Mom.mean_lin are such that Mom.mean_cost + Mom.mean_lin*x0  = Mom.mean        

% Views = struct('N_MeanViews','N_CovViews','dimension','monitoring_time','v_mu','v_sig','mu_view','sig2_view')
% where:
% N_MeanViews = Number of views on the expectations              [scalar]
% N_CovViews = Number of views on the covariance matrix          [scalar]
% dimension = labels of the risk drivers                         [t_*n_ x 1]
% monitoring_time = monitoring times                             [t_*n_ x 1]
% v_mu = matrix that qualifies the views on expectations         [N_MeanViews x t_*n_]
% v_sig = matrix that qualifies the views on covariance          [N_CovViews x t_*n_]
% mu_view = extent of the views on expectation                   [N_MeanViews x 1]
% sig2_view = extent of the views on the covariance              [N_CovViews x 1]

% OUTPUT
% Posterior = struct('monitoring_time','dimension','mean', 'cov','mean_cost','mean_lin');
% where:
% Posterior.monitoring_time = monitoring times                        [t_*n_ x 1]
% Posterior.dimension = labels of the risk drivers                    [t_*n_ x 1]    
% Posterior.cov = posterior covariance matrice of X_{t1,t2...t_}      [t_*n_ x t_*n_]
% Posterior.mean = posterior vector of the means of X_{t1,t2...t_}    [t_*n_ x 1]
% Posterior.mean_cost                                                 [t_*n_ x 1]
% Posterior.mean_lin                                                  [t_*n_ x 1]
% Posterior.mean_cost and Posterior.mean_lin are such that Posterior.mean_cost + Posterior.mean_lin*x0  = Posterior.mean        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = unique(Mom.dimension);
t = unique(Mom.monitoring_time);
n_ = length(n);
t_ = length(t);
[T,N] = meshgrid(t,n);
Posterior = struct('monitoring_time', T(:), 'dimension',N(:),'mean', NaN(t_*n_,1), 'cov',NaN(n_*t_,n_*t_),'mean_cost',NaN(t_*n_,1),'mean_lin',NaN(t_*n_,n_));

S2 = Mom.cov;
mu = Mom.mean;
v_mu = Views.v_mu;
v_sig = Views.v_sig;
mu_view = Views.mu_view;
sig2_view = Views.sig2_view;

if isnan(Views.v_mu) == 1
    Posterior.mean = mu;
else
    v_mu_dag = (S2*v_mu')/(v_mu*S2*v_mu');
    P_orth = v_mu_dag*v_mu;
    P = eye(size(P_orth)) - P_orth;    
    Posterior.mean = P*mu+P_orth*v_mu_dag*mu_view;
    Posterior.mean_lin = P*Mom.mean_lin;
    Posterior.mean_cost = P*Mom.mean_cost + P_orth*v_mu_dag*mu_view;
end
if isnan(Views.v_sig) == 1
    Posterior.cov = S2;
else    
    v_sig_dag = (S2*v_sig')/(v_sig*S2*v_sig');
    P_orth = v_sig_dag*v_sig;
    P = eye(size(P_orth)) - P_orth;       
    Posterior.cov = P*S2*P'+P_orth*(v_sig_dag*sig2_view*v_sig_dag')*P_orth';
end

