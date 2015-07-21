function q_t = QuadraticMat_Vc(lambda,gamma,eta,sig2,c2,tau_,n_,i_invest)

% This function computes the matrix q_t of the problem to solve when using CALCULUS of VARIATION:
% argmin_b (b' q_t b - b'l_t)

% INPUT
% lambda = discounting parameter                                        [scalar]
% gamma = risk aversion parameter                                       [scalar]
% eta = overall weight of the market impact of transactions             [scalar]
% sig2 = covariance matrix of the process of the risk drivers           [n_*t_ x n_*t_]
% c2 = market impact matrix                                             [k_ x k_]
% tau_ = effective number of future time steps considered               [scalar]
% n_ = number of risk drivers                                           [scalar]
% i_invest = labels of the investible risk drivers                      [k_ x 1]
% where: 
% t_ = number of monitoring times at which sig2 is computed 
% k_ = number of investible risk drivers

% OUTPUT
% q_t = matrix                                                          [k_*tau_ x k_*tau_]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin  < 8 || isempty(i_invest)
    i_invest = 1:1:length(c2);
end

k_ = length(i_invest);%number of investible risk drivers
for t = 1:tau_-1
    sig2_t = ExtractBlockMtx_fn(sig2,t,t,n_,n_);
    sig2_t = sig2_t(i_invest,i_invest);
    sig2_t1 = ExtractBlockMtx_fn(sig2,t+1,t+1,n_,n_);
    sig2_t1 = sig2_t1(i_invest,i_invest);
    sig2_tt1 = ExtractBlockMtx_fn(sig2,t,t+1,n_,n_);
    sig2_tt1 = sig2_tt1(i_invest,i_invest);
    sig2_t = sig2_t + sig2_t1 -2*sig2_tt1; 
    q_t((t-1)*k_+1:t*k_,(t-1)*k_+1:t*k_) = exp(-lambda*(t-1))*(-gamma*0.5*sig2_t-eta*0.5*c2*(1+exp(-lambda)));   
    q_t((t-1)*k_+1:t*k_,t*k_+1:(t+1)*k_) = exp(-lambda*t)*eta*0.5*c2;       
    q_t(t*k_+1:(t+1)*k_,(t-1)*k_+1:t*k_) = exp(-lambda*t)*eta*0.5*c2;       
end

t = tau_;
sig2_t = ExtractBlockMtx_fn(sig2,t,t,n_,n_);
sig2_t = sig2_t(i_invest,i_invest);
sig2_t1 = ExtractBlockMtx_fn(sig2,t+1,t+1,n_,n_);
sig2_t1 = sig2_t1(i_invest,i_invest);
sig2_tt1 = ExtractBlockMtx_fn(sig2,t,t+1,n_,n_);
sig2_tt1 = sig2_tt1(i_invest,i_invest);
sig2_t = sig2_t + sig2_t1 -2*sig2_tt1; 
q_t((t-1)*k_+1:t*k_,(t-1)*k_+1:t*k_) = exp(-lambda*(t-1))*(-gamma*0.5*sig2_t-eta*0.5*c2);
q_t = (q_t + q_t')/2;



