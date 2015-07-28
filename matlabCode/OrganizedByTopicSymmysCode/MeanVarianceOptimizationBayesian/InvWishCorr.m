function av_cor = InvWishCorr(S,Nu,NumCases)

N=size(S,1);
%Mult=2*N/(N-1);
invSNu=inv(S)/Nu;

av_cor=[];
for n=1:NumCases
  % generate W ~ W(Nu,inv(S)/Nu);
  X = mvnrnd(zeros(N,1),invSNu,Nu);
  W=0*S;
  for n=1:Nu
    W=W+X(n,:)'*X(n,:);
  end
  % compute average correlation (modulo constant) of inv(W) 
  [a,c]=cov2corr(inv(W));
  av_cor=[av_cor sum(sum(c))];
end
av_cor=(av_cor-N)/(N*(N-1));

