function [mu, ATMF, BTMF, CTMF, Hhat] = FlexM(returns,demean,eps,df)
% PURPOSE: estimate FlexM multivariate GARCH model
% ------------------------------------------------------------------
% Where: returns is a matrix of asset returns of dimension T x N;
%                so rows must correspond to time and columns to assets
%        demean specifies whether returns should be demeaned (if demean = 1)
%               or not (otherwise) to estimate the model; default value is 1
%        eps is used in enforcing a_ii + b_ii <= 1 - eps;
%               the default value is zero
%        df is the degree of freedom for the t-distribution;
%               the default value is 500 to make it, basically, normal
% ------------------------------------------------------------------
% RETURNS: ATMF = coefficient matrix A-tilde (in the notation of the paper)
%          BTMF = coefficient matrix B-tilde (in the notation of the paper)
%          CTMF = coefficient matrix C-tilde (in the notation of the paper)
%          Hhat = forecasted conditional covariance matrix
% ------------------------------------------------------------------
% NOTES: 
% ------------------------------------------------------------------

% written by:
% Olivier Ledoit and Michael Wolf
% CREATED  04/18/02
% UPDATED 

% default value for df is 500
if nargin < 4
  df = 500;
end

% default value for eps is 0
if nargin < 3
  eps = 0;
end

% default value for demean is 1
if nargin < 2
   demean = 1;
end

if (eps < 0) 
  error('eps must be a (small) positive number')
end


% Initialization
[T,N]=size(returns);
if (1 == demean)
   mu = mean(returns);
   returns = returns-mu(ones(1,T),:);
end
S = returns'*returns/(T-1);
x = returns';

A=zeros(N,N);
B=zeros(N,N);
C=zeros(N,N);

% Rescale Data
scale=sqrt(mean(x'.^2)');
x=x./scale(:,ones(T,1));

% Estimation of On-Diagonal Elements
h=zeros(N,T);
for i=1:N
   % Likelihood Maximization
   q=garch1f4(x(i,:)',eps,df);
	A(i,i)=q(2);
	B(i,i)=q(3);
	C(i,i)=q(1);
	h(i,:)=filter([0 q(2)],[1 -q(3)],x(i,:).^2*(df-2)/df,mean(x(i,:).^2)*(df-2)/df)...
	       +filter([0 q(1)],[1 -q(3)],ones(1,T));
end

% First-step Estimation of Off-Diagonal Elements
for i=1:N
	for j=i+1:N
      % Likelihood Maximization
		theta=garch2f8(x(i,:).*x(j,:),C(i,i),A(i,i),B(i,i),x(i,:).^2,h(i,:),C(j,j),A(j,j),B(j,j),x(j,:).^2,h(j,:),df);
      A(i,j)=theta(2);
		B(i,j)=theta(3);
		C(i,j)=theta(1);
		A(j,i)=A(i,j);
		B(j,i)=B(i,j);
		C(j,i)=C(i,j);
   end
end

% Transformation of Coefficient Matrices
ATMF=minfro(A);
BTMF=minfro(B);
CTMF=minfro(C./(1-B)).*(1-BTMF);

% Rescale
C=C.*(scale*scale');
CTMF=CTMF.*(scale*scale');

% Forecast of Conditional Covariance Matrix
Hhat = zeros(N,N);
for i = 1:N
   for j = 1:N
      hSeries = filter([0 ATMF(i,j)], [1 -BTMF(i,j)], returns(:,i).*returns(:,j), S(i,j))...
         + filter([0 CTMF(i,j)], [1 -BTMF(i,j)], ones(T,1));
      Hhat(i,j) = hSeries(T);
   end
end

mu = mu';



