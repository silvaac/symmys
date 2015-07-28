function [b_prior b_post] = BellmanEq_CS1(eta, gamma, lambda, tau, theta, mu, sig2, c2, b_legacy, x, t_view, view)

% This function solves the Bellman Equation for the case study 1. 
% In Case Study 1 there is only one risk driver (the 10 year rate) and only 
% one view. The view is that the expected value of the 10-year rate will be 
% the actual value minus 50 basis points at t^view = 1 year from the current time.
% The solution is analytical in the prior case.
% The solution is found recursivelly in the posterior case starting from
% t_view and going back to time 0 with a time step = tau.

% INPUT
% eta = overall weight of the market impact of transactions     [scalar]
% gamma = risk aversion parameter                               [scalar]
% lambda = discounting parameter                                [scalar]
% tau = trading interval                                        [scalar]
% theta = transition matrix of the MVOU process                 [n_ x n_]
% mu = drift vector of the MVOU process                         [n_ x 1]
% sig2 = covariance parameters of the MVOU process              [n_ x n_]
% c2 = matrix of the market impact                              [n_ x n_]
% b_legacy = legacy portfolio exposure at time 0                [n_ x 1]
% x = path of the risk drivers (with time step = tau)           [t_ x n_]
% t_view = times of the views                                   [1 x N_MeanViews] 
% view = views on the risk drivers                              [1 x N_MeanViews]
% In case study 1: n_ = 1; N_meanViews = 1

% OUTPUT
% b_prior = optimal prior exposure                              [t_x n_ matrix]                    
% b_post = optimal posterior exposure                           [t_x n_ matrix] 


t_ = length(x);
n_ = length(theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the prior at time 0
Prior0 = MVOU_Prior([0 tau], x(1), theta, sig2, mu);
% first period covariance matrix
sig2_1 = Prior0.cov(n_+1:2*n_,n_+1:2*n_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coefficients of the Bellman equation according to the prior 

alpha_prior = Prior0.mean_cost(n_+1:2*n_);
beta_prior = Prior0.mean_lin(n_+1:2*n_,1:n_)-eye(n_);
HATsig2 = exp(lambda)*(eta*c2)^(-1/2)*gamma*sig2_1*(eta*c2)^(-1/2);
HATpsi_bb = (0.25*(HATsig2+eye(n_)*(exp(lambda)-1))^2 + HATsig2)^(1/2) - 0.5*(HATsig2+eye(n_)*(exp(lambda)-1));
psi_bb_prior = (eta*c2)^(1/2)*HATpsi_bb*(eta*c2)^(1/2);
q_prior = gamma*sig2_1 + eta*c2 + exp(-lambda)*psi_bb_prior;
tmp = (eta*c2*(q_prior\beta_prior));
psi_bx_prior = (eye(n_^2)-exp(-lambda)*kron(beta_prior'+eye(n_),eta*(c2/q_prior)))\tmp(:);
psi_bx_prior = reshape(psi_bx_prior,n_,n_);
psi_b_prior = (q_prior/(eta*c2) - exp(-lambda)*eye(n_))\(eye(n_)+exp(-lambda)*psi_bx_prior)*alpha_prior;

% %Alternatively, if and only if c2 = sig2_1
% a_ = (sqrt(4*gamma*eta*exp(-lambda)+(gamma+(1-exp(-lambda))*eta)^2) - ((1-exp(-lambda))*eta+gamma))/(2*exp(-lambda));
% psi_bb_prior = a_*sig2_1;
% psi_bx_prior = eta*(beta_prior)*inv((gamma+eta+exp(-lambda)*a_)*eye(n_)-eta*exp(-lambda)*(beta_prior+eye(n_)));
% psi_b_prior = eta*(alpha_prior + exp(-lambda)*psi_bx_prior*alpha_prior)/(gamma+eta+exp(-lambda)*a_-exp(-lambda)*eta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coefficients of the Bellman equation according to the posterior distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inizialize the variables
Hor = ceil(max(t_view)/tau);
psi_t_bb = zeros(n_,n_,Hor);
q_t = zeros(n_,n_,Hor);
psi_t_bx = zeros(n_,n_,Hor);
psi_t_b = zeros(n_,Hor);
alpha_t = zeros(n_,Hor);
beta_t = zeros(n_,n_,Hor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set the boundary conditions asintotically. After the last view, the
%solution is equal to the prior
psi_t_bb(1:n_,1:n_,Hor) = psi_bb_prior; 
q_t(1:n_,1:n_,Hor) = q_prior;
psi_t_bx(1:n_,1:n_,Hor) = psi_bx_prior;
psi_t_b(1:n_,Hor) = psi_b_prior;
alpha_t(1:n_,Hor) = alpha_prior;
beta_t(1:n_,1:n_,Hor) = beta_prior;

for k = Hor-1:-1:1    
    q = squeeze(q_t(1:n_,1:n_,k+1));
    alpha = squeeze(alpha_t(1:n_,k+1));
    beta = squeeze(beta_t(1:n_,1:n_,k+1));
    psi_bx = squeeze(psi_t_bx(1:n_,1:n_,k+1));
    psi_b = squeeze(psi_t_b(1:n_,k+1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update of the coefficients of the Bellman equation
    psi_t_b(1:n_,k) = eta*c2*(q\(alpha + exp(-lambda)*psi_bx*alpha+exp(-lambda)*psi_b));
    psi_t_bx(1:n_,1:n_,k) = eta*c2*(q\(beta + exp(-lambda)*psi_bx*(beta+eye(n_)))); 
    psi_t_bb(1:n_,1:n_,k) = -eta^2*c2*(q\c2) + eta*c2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %set the monitoring times of interest
    t = [0 tau t_view - k*tau];
    if (t(end)-t(end-1))<t(end)*10^-10
        t = t(1:end-1);
    end
    T_ = length(t);

    %compute the prior
    Prior = MVOU_Prior(t, 0, theta, sig2, mu);    
    
    %set the views
    [T,N] = meshgrid(t,1:n_);
    N_Meanviews = 1;%Number of views on expectations
    v_mu_tmp = zeros(N_Meanviews,n_,T_);
    v_mu_tmp(1,1,T_) = 1;
    mu_view(1,1) = view(1); 
    v_mu = v_mu_tmp(:,:);
    views = struct('N_Meanviews', N_Meanviews,'N_Covviews',[],'dimension',N(:),'monitoring_time',T(:),...
        'v_mu',v_mu,'v_sig',NaN,'mu_view',mu_view,'sig2_view',[]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update the Posterior moments of the process
    Posterior = MVOU_Posterior(Prior, views);    
    alpha_t(1:n_,k) = Posterior.mean_cost(n_+1:2*n_);
    beta_t(1:n_,1:n_,k) = Posterior.mean_lin(n_+1:2*n_,1:n_)-eye(n_);    
    sig2_1 = Posterior.cov(n_+1:2*n_,n_+1:2*n_);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update matrix q_t
    psi_bb = squeeze(psi_t_bb(1:n_,1:n_,k));
    q_t(1:n_,1:n_,k) = gamma*sig2_1 + eta*c2 + exp(-lambda)*psi_bb;    
end

b_prior = NaN(t_,1);
b_post = NaN(t_,1);

%Reconstructing the optimal exposure on the simulated path x
for t = 1:t_
    if t == 1
        b_legacy_prior = b_legacy;
        b_legacy_post = b_legacy;
    else
        b_legacy_prior = b_prior(t-1);
        b_legacy_post = b_post(t-1);        
    end
    % prior exposure 
    b_prior(t) = q_prior\(alpha_prior+(beta_prior)*x(t)+eta*c2*b_legacy_prior...
        +exp(-lambda)*psi_bx_prior*(alpha_prior+(beta_prior+eye(n_))*x(t))+exp(-lambda)*psi_b_prior);        
    %posterior exposure
    q = squeeze(q_t(1:n_,1:n_,t));
    alpha = squeeze(alpha_t(1:n_,t));
    beta = squeeze(beta_t(1:n_,1:n_,t));
    psi_bx = squeeze(psi_t_bx(1:n_,1:n_,t));
    psi_b = squeeze(psi_t_b(1:n_,t)); 
    l = alpha + beta*x(t)+eta*c2*b_legacy_post+exp(-lambda)*psi_bx*(alpha+(beta+eye(n_))*x(t))+exp(-lambda)*psi_b; 
    b_post(t) = q\l;
end





