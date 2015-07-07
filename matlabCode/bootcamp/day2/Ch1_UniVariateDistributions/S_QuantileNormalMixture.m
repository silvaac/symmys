function Q=QuantileNormalMixture(c,ps,mus,sigs)
% this function computes the quantile for a mixture of normals 
% the confidence level p can be a Jx1 vector. If this vector is uniformly
% distributed on [0,1] the ensuing random variable Q is distributed as the mixture of normals
% for more information, see "Risk and Asset Allocation"-Springer (2005), by A. Meucci


% compute first moment
m=ps'*mus; 
% compute second moment
E2=ps'*(mus.^2+sigs.^2);
s=sqrt(E2-m^2);

% compute cdf on suitable range
Z = m + 6*s*[-1: .001 : 1];

F = 0*Z;
for n=1:length(mus)
    F=F+ps(n)*normcdf(Z,mus(n),sigs(n));
end

% compute quantile by interpolation
Q = interp1(F,Z,c,'linear','extrap');

    