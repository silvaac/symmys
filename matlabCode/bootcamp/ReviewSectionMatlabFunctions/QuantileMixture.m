function Q=QuantileMixture(p,a,m_Y,s_Y,m_Z,s_Z)
% this function computes the quantile of a mixture distirbution by linear interpolation/extrapolation of the cdf
% the confidence level p can be vector. If this vector is uniformly distributed on [0,1] 
% the sample Q is distributed as the mixture
% see Section 3.4.5 of "Risk and Asset Allocation" - Springer (2005), by A. Meucci

% compute first moment
m=a*m_Y+(1-a)*exp(m_Z+.5*s_Z*s_Z); 
% compute second moment
Ex2=a*(m_Y^2+s_Y^2)+(1-a)*exp(2*m_Z+2*s_Z*s_Z);
s=sqrt(Ex2-m*m);

% compute cdf on suitable range
X = m + 6*s*[-1: .001 : 1];
F = a*normcdf(X,m_Y,s_Y)+(1-a)*logncdf(X,m_Z,s_Z);
[F, ind] = unique(F);
X=X(ind);

% compute quantile by interpolation
Q = interp1(F,X,p,'linear','extrap');

    