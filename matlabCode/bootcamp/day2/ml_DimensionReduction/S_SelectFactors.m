% this script selects the best stocks to replicate an index. The same routine
% works more in general to select the best explicit factors from a large pool of potential candidates.
% see Section 3.4.5 of "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; close all; clear;

% load database of Stock prices P, S&P500 M, and dates D 
load DatabaseStocks

% compute returns (single stocks and S&P)
R=P(2:end,:)./P(1:end-1,:);
R_M=M(2:end)./M(1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-select factors by picking one stock per correlation cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=cov([R R_M]);
s=sqrt(diag(S));
C=diag(1./s)*S*diag(1./s);

% computes distance ~ according to square correlation
Pdist = [];
for n=1:size(R,2)-1
    Pdist = [Pdist 1-(C(n,n+1:end-1)).^2];
end
Z = linkage(Pdist ,'complete'); 
T = cluster(Z,'cutoff',1.15);

% in each cluster keep the stock that is most correlated with the market
Keep=[];
for n=1:max(T)
    ClusterIndex=find(T==n);
    CorrClusterMarket=C(ClusterIndex,end);
    [dummy,IndexMaxCorr]=max(CorrClusterMarket.^2);
    Keep=[Keep
        ClusterIndex(IndexMaxCorr)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine the best pool of stocks for a given size K with the routine in
% Section 3.4.5 of "Risk and Asset Allocation" - Springer (2005), by A. Meucci
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define dependent variable(s) X (this can be multivariate) in X=BF+U
X = R_M;
L = size(X,2);   % dimension of invariants X in X=BF+U

% define replicating factors in X=BF+U (see clustering routine above)
F=R(:,Keep);
N=size(F,2);   % number of factors F in X=BF+U

Y=[X F];
cov_Y=cov(Y);  % compute joint covariance of invariants and explanatory factors

Cov_XX=cov_Y(1:L,1:L);         % extract sub-covariance of invariants
Cov_FF=cov_Y(L+1:end,L+1:end); % extract sub-covariance of factors
Cov_XF=cov_Y(1:L,L+1:end);     % extract cross-covariances of invariants and factors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here starts the routine
K=N;
I_K=[1:N];  % initialize pool of all factors
v_K_star=[];
I_K_star=[];
NumFactors=[];

while K>1
    CyclesToGo=K
    v_k=[];
    for k=1:K
        I_k=I_K; I_k(k)=[];                       % drop the k-th element from the pool I_K
        v_k(k)=R2(I_k,Cov_XX,Cov_XF,Cov_FF);      % evaluate target (R^2 in this case)
    end
    [v_k_star,k_star]=min(v_k);      % find best element
    I_K(k_star)=[];                  % get rid of best element
    
    K=length(I_K);
    
    v_K_star=[v_K_star v_k_star];    % record the optimal target at each iteration (R^2 in this case)
    NumFactors=[NumFactors K];
    I_K_star(K).Record=I_K;        % record the optimal pool of size K
                                   % to convert to stock index use Keep(I_K_star(K).Record)
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the optimal target (R^2) as function of the iteration
figure
plot(NumFactors,v_K_star)    
xlabel('number of stocks')
ylabel('r-square')
grid on