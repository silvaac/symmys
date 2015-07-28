function XXX=minfro(A);

% function XXX=minfro(A);
%
% Copyright Ilya Sharapov (1997)
% Input:  Matrix A is an indefinite symmetric matrix with 
%         non-negative diagonal elements
% Output: Matrix XXX is the positive semi-definite matrix
%         with same diagonal elements as A that is closest
%         to A according to the Frobenius norm

if any(diag(A)<0)
	error('Diagonal Elements Must Be Non-Negative!')
elseif any(any(A~=A'))
	error('Matrix Must Be Symmetric!')
elseif all(eig(A)>=0)
	XXX=A;
else

rho=0.75;		% IMPORTANT parameter: 
			% if things go wrong make rho 
			% bigger and wait longer 

tol=3e-6; 		% tolerance 
maxj=10;		% max number of iterations 

[n,nn]=size(A);
M=diag(diag(A));	% initialize with diagonal 
[n,nn]=size(A);
oldnorm=norm(M-A,'fro'); 
oldnormj=oldnorm;
normj(1)=oldnorm;

j=1;
incmax=1e32;		% just to enter the loop 
while ((j<maxj)&(incmax>tol))
   incmax=0; 
%fprintf('\n Iteration Subiteration Objective Increment Infeasability \n'); 
%fprintf('  --------------------------------------------------------- \n');
   for i=1:n;
      a=[A(1:i-1,i);A(i+1:n,i)];
      m=[M(1:i-1,i);M(i+1:n,i)];
      aii=A(i,i);
      b=a-rho*m;
      x=newton(M,i,b,m,aii,n,rho);	% Newton's step 
      P=sparse(eye(n));
      P(i,1:n)=x'; 
      Mtest=P*M*P';			% update 
      M=Mtest;
      inc=oldnorm-norm(M-A,'fro');
      oldnorm=norm(M-A,'fro');
      if (inc > incmax) 		% find maximal increment 
         incmax=inc;			% over iteration
      end
      infeas=max(abs(diag(M)-diag(A)));
%      fprintf('%8i %8i %16.4e %8.4e %8.4e  \n', j, i, oldnorm, inc, infeas);
   end
   normj(j+1)=oldnorm;
   incj(j)=oldnormj-oldnorm;
   oldnormj=oldnorm;
%   fprintf('\n Iteration %4i completed \n', j) ;
%   fprintf(' Maximal increment for subiteration %8.4e \n', incmax);
%   fprintf(' Increment over iteration %8.4e \n\n', incj(j));
   j=j+1;
end

XXX=M;

end