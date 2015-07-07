function F_U = TCopulaCDF(u,nu,Mu,Sigma)

N=length(u);
[s,C] = cov2corr(Sigma);

x=Mu+s'.*tinv(u,nu);

F_U = real(mvtcdf(x-Mu./s',C,nu));


