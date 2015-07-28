clear all
% Case Study 2
% In Case Study 2 we consider two risk drivers, the 10 year rate and the TIP spread, and two 
% non-synchronous views on them. The view on the rate is that its expected 
% value will be the actual value minus 50 basis points
% at t_viewX = 1 year from the current time (as in Case Study 1). The view
% on the TIP spread is that its expected value will be the actual value
% plus 50 basis points at t_viewTIP = 0.75 years.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the Calibrated Parameters
load 'CalibratedParameters.mat'
mu_LT = theta\mu;
% 
%load the simulated path for the 10y Government rate (this is equal to its
%shadow rate) and the TIP spread
load 'path_daily_rate&TIP1'
x0 = [r_path(1); TIP_path(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set the trading
tau = 1/252;                        % trading frequency (daily)
T_Hor = 0.5;%years                  % effective future portfolio horizon   
tau_ = length([0:tau:T_Hor]);       % effective number of future tradings at any point in time
n_ = 2;                             % number of risk drivers 
i_invest = 1;                       % the first of the risk drivers (the 10y rate) is the investible one
k_ = length(i_invest);              % number of investible risk drivers
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set the view at time 0
Index_t_viewX = ceil(1/tau);    
t_viewX = tau*Index_t_viewX;%years                          % time of the view on the rate 
mu_X = x0(1) - 0.5/100;                                     % view on the 10y rate   
Index_t_viewTIP = ceil(0.75/tau);
t_viewTIP = tau*Index_t_viewTIP;%years                      % time of the view on the TIP spread
mu_TIP = x0(2)+0.5/100;                                     % view on the TIP spread
t = sort(unique([0:tau:T_Hor t_viewX t_viewTIP]'));         % monitoring times
t_ = length(t);                                             % number of monitoring times    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set the parameters for optimization
gamma = 10^-2;                                              % risk aversion parameter    
eta = 0.5;                                                  % weight of the market impact of transaction
lambda = log(2)/20;                                         % discount (half life 20*tau)
b_legacy = 0;                                           % legacy portfolio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the prior at time 0
Prior0 = MVOU_Prior(t, x0, theta, sig2, mu);
%matrix of market impact
c2 = Prior0.cov(n_+1:2*n_,n_+1:2*n_);
c2 = c2(i_invest,i_invest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prior and Posterior Optimal exposure with Market Impcat of transaction.
% SOLVING THE BELLMAN EQUATION

t_view0 = [t_viewX t_viewTIP];
view = [mu_X mu_TIP];
i_view = [1 2];
omega = [1,0];
[b_MI_Bellman_prior b_MI_Bellman_post] = BellmanEq_CS2(eta, gamma, lambda, tau, theta, mu, sig2, c2, b_legacy, X_path, t_view0,view,i_view,omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULUS OF VARIATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize variables

%prior
b_NoMI_prior = NaN(tau_-1,k_);      %optimal 1-period myopic prior solution 
b_MI_Vc_prior = NaN(tau_-1,k_);     %optimal prior solution with MI (Calculus of Variation)
b_legacy_prior = b_legacy;
b_legacy_prior_LR = b_legacy;
b_legacy_prior_xt = b_legacy;

%posterior 
b_NoMI_post = NaN(tau_-1,k_);           %optimal 1-period myopic posterior solution
b_NoMI_LongTermX = NaN(tau_-1,k_);      %Decomposition of the optimal solution in absence of market impact
b_NoMI_viewMeanX = NaN(tau_-1,k_);      %Decomposition of the optimal solution in absence of market impact
b_NoMI_LongTermTIP = NaN(tau_-1,k_);    %Decomposition of the optimal solution in absence of market impact
b_NoMI_viewMeanTIP = NaN(tau_-1,k_);    %Decomposition of the optimal solution in absence of market impact
b_MI_Vc_post = NaN(tau_-1,k_);          %optimal posterior solution with MI (Calculus of Variation)
b_legacy_post = b_legacy;
b_legacy_post_LR = b_legacy;
b_legacy_post_xt = b_legacy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_tmp = 1;
for i = 1:tau_-1
    i
    if  (t_viewTIP - tau*(i-1) > T_Hor)
        t_roll = [[0:tau:T_Hor] t_viewTIP-tau*(i-1) t_viewX-tau*(i-1)]';
        Index_t_viewTIP = t_-1;
    else
        t_roll = [[0:tau:T_Hor]  t_viewX - tau*(i-1)]';
        Index_t_viewTIP = length(0:tau:T_Hor)-(i_tmp-1);
        i_tmp = i_tmp+1;
    end
    t_ = length(t_roll);
    Index_t_viewX = t_;
    t_viewTIP_roll = t_roll(Index_t_viewTIP);
    t_viewX_roll = t_roll(Index_t_viewX);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                   The Prior
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %compute the prior
    Prior = MVOU_Prior(t_roll, X_path(i,:)', theta, sig2, mu);    
    Mean_prior = reshape(Prior.mean(1:t_*n_),n_,t_)';     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% No Market impact myopic 1 period solution
    ER_prior = Mean_prior(2,i_invest)-Mean_prior(1,i_invest);       %expected return of the 10y-rate
    sig2_1 = Prior.cov(n_+1:2*n_,n_+1:2*n_);
    sig2_1 = sig2_1(i_invest,i_invest);                             % variance of the 10y-rate  
    b_NoMI_prior(i,i_invest) = (gamma*sig2_1)\(ER_prior(1:k_)');    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Market impact: calculus of variation
    ER = diff(Mean_prior(1:end-2,i_invest));
    l_t = - exp(-lambda*[0:1:length(ER)-1]').*ER;
    l_t(1) = l_t(1) - eta*c2*b_legacy_prior;
    q_t = QuadraticMat_Vc(lambda,gamma,eta,Prior.cov,c2,length(ER),n_,i_invest);
    b_MI_Vc_prior_tmp = (2*q_t)\l_t;
    b_MI_Vc_prior(i,1:k_) = b_MI_Vc_prior_tmp(1:k_);
    b_legacy_prior = b_MI_Vc_prior(i,1:k_);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%The Posterior
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %set the rolling view 
    [T,N] = meshgrid(t_roll,1:n_);
    clear v_mu v_tmp mu_view;
    N_Meanviews = 2;%Number of views on expectations
    v_tmp = zeros(N_Meanviews,n_,t_);
    %first view (on the rate)
    v_tmp(1,1,Index_t_viewX) = 1; %view on the rate
    mu_view(1,1) = mu_X; 
    %second view (on TIP spread)
    v_tmp(2,2,Index_t_viewTIP) = 1; %view on inflation
    mu_view(2,1) = mu_TIP; 
    v_mu = v_tmp(:,:);
    views = struct('N_Meanviews', N_Meanviews,'N_Covviews',[],'dimension',N(:),'monitoring_time',T(:),...
        'v_mu',v_mu,'v_sig',NaN,'mu_view',mu_view,'sig2_view',[]);                    
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the posterior moments
    Posterior = MVOU_Posterior(Prior, views);   
    Mean_post = reshape(Posterior.mean(1:t_*n_),n_,t_)';     

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% No Market impact myopic 1 period solution
    ER_post = Mean_post(2,i_invest)-Mean_post(1,i_invest);          %expected return of the 10y-rate
    sig2_1 = Prior.cov(n_+1:2*n_,n_+1:2*n_);
    sig2_1 = sig2_1(i_invest,i_invest);                             % variance of the 10y-rate  
    b_NoMI_post(i,i_invest) = (gamma*sig2_1)\(ER_post(1:k_)');    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Market impact: calculus of variation
    ER = diff(Mean_post(1:end-2,i_invest));
    l_t = - exp(-lambda*[0:1:length(ER)-1]').*ER;
    l_t(1) = l_t(1) - eta*c2*b_legacy_post;
    q_t = QuadraticMat_Vc(lambda,gamma,eta,Prior.cov,c2,length(ER),n_,i_invest);
    b_MI_Vc_post_tmp = (2*q_t)\l_t;
    b_MI_Vc_post(i,1:k_) = b_MI_Vc_post_tmp(1:k_);
    b_legacy_post = b_MI_Vc_post(i,1:k_);
         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% LongTerm and viewMean contributions of the optimal solution with views in absence of market impact of transactions 
        
    theta1 = theta(1,1);
    theta2 = theta(2,2);

    %extract the covariance matrix relative to time t, t+1, t_viewTIP_roll and t_viewX_roll 
    idx = [1,2,3,4, n_*Index_t_viewTIP-1, n_*Index_t_viewTIP, t_*n_-1, t_*n_];
    SubSig2 = Posterior.cov(idx,idx);   
    delta = (SubSig2(3,7)*SubSig2(6,6)-SubSig2(3,6)*SubSig2(6,7))/(SubSig2(6,6)*SubSig2(7,7)-SubSig2(6,7)^2);
    rho = (SubSig2(3,6)*SubSig2(7,7)-SubSig2(3,7)*SubSig2(6,7))/(SubSig2(6,6)*SubSig2(7,7)-SubSig2(6,7)^2);
    
    b_NoMI_LongTermX(i,1:k_) = 2*theta1/(gamma*sig2(1,1))/(1+exp(-theta1*tau))*(1-delta*(1-exp(-theta1*t_viewX_roll))/(1-exp(-theta1*tau)))*(mu_LT(1)-X_path(i,1)); 
    b_NoMI_LongTermTIP(i,1:k_) = -2*rho*theta1/(gamma*sig2(1,1))/(1-exp(-2*theta1*tau))*(1-exp(-theta2*t_viewTIP_roll))*(mu_LT(2)-X_path(i,2));
    
    b_NoMI_viewMeanX(i,1:k_) = 2*delta*theta1/(gamma*sig2(1,1))/(1-exp(-2*theta1*tau))*(mu_view(1) - X_path(i,1));%(view on X)
    b_NoMI_viewMeanTIP(i,1:k_) = 2*rho*theta1/(gamma*sig2(1,1))/(1-exp(-2*theta1*tau))*(mu_view(2) - X_path(i,2));%(view on TIP)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save ('data\CaseStudy2.mat','X_path','r_path','TIP_path','t','mu_LT','mu_view','b_NoMI_prior','b_NoMI_post',...
%     'b_MI_Vc_prior','b_MI_Vc_post','b_MI_Bellman_prior','b_MI_Bellman_post',...
%     'b_NoMI_LongTermX','b_NoMI_LongTermTIP','b_NoMI_viewMeanX','b_NoMI_viewMeanTIP','tau_');

