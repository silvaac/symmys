function xx=newton(M,i,b,m,aii,n,rho);

% function xx=newton(M,i,b,m,aii,n,rho);
%
% Copyright Ilya Sharapov (1997)
% Subroutine called interbally by minfro.m

maxit=40;
eps=1e-9;               % small correction 
			% to handle singularity 
l=0.0;
MM=[M(1:i-1,1:i-1) M(1:i-1,i+1:n); M(i+1:n,1:i-1) M(i+1:n,i+1:n)]+eps*eye(n-1);

j=1;
while (j<maxit)           % loop         
   IM=inv(MM*MM+l*MM);
   x=IM*(MM*b-l*rho*m);
   f=rho*rho*aii+2*rho*x'*m+x'*MM*x-aii;
   if (abs(f)<1e-7)       %
      break;              % stopping
   end                    %
   dfdl=-2*(rho*m+MM*x)'*IM*(rho*m+MM*x);
   l=l-f/dfdl;              % Newton's step
   j=j+1;
end
if (abs(f)<1e-7)			% converged 
   xx=[x(1:i-1); rho; x(i:n-1)];
else					% didn't converge
   xx=zeros(n,1);			%
   xx(i)=1;				% this vector will not 
   xx;					% change M 
end


