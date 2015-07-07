function r2=R2(I_k,Cov_XX,Cov_XF,Cov_FF)

Cov_FF_k=Cov_FF(I_k,I_k);
Cov_XF_k=Cov_XF(:,I_k);

Cov_U=Cov_XX-Cov_XF_k*(Cov_FF_k\(Cov_XF_k'));

r2=1-trace(Cov_U)/trace(Cov_XX);