function q=garch2f8(y,c1,a1,b1,y1,h1,c2,a2,b2,y2,h2,df)

% function q=garch2f7(y,c1,a1,b1,y1,h1,c2,a2,b2,y2,h2)
%
% Off-diagonal parameter estimation in bivariate GARCH(1,1)
% when diagonal parameters are given.
%
% Uses a conditional t-distribution with fixed degrees of freedom
%
% Steepest Ascent on boundary
% Hessian off boundary
% No grid search

% Parameters
gold=(1+sqrt(5))/2;	% step size increment 
tol1=1e-7; 		% for termination criterion
tol2=1e-7;		% for closeness to boundary
big=2;			% for making the hessian negative definite
maxiter=50;		% maximum number of iterations
n=30;			% number of points on the grid

% Prepare
t=length(y);
y1=y1(:);
y2=y2(:);
y=y(:);
s=mean(y);
s1=mean(y1);
s2=mean(y2);
h1=h1(:);
h2=h2(:);

% Bounds
low=[-sqrt(c1*c2) 0 0]+tol2;
high=[sqrt(c1*c2) sqrt(a1*a2) sqrt(b1*b2)]-tol2;

% Starting Point
a0=0.9*sqrt(a1*a2);		
b0=0.9*sqrt(b1*b2);		
c0=mean(y)*(1-a0-b0)*(df-2)/df;
c0=sign(c0)*min(abs(c0),0.9*sqrt(c1*c2));
      
      % Initialize optimization
a=[c0 a0 b0];
best=0;
da=0;
term=1;
negdef=0;
iter=0;

% Begin optimization loop
while iter<maxiter;
iter=iter+1;

% New parameter
olda=a;
a=a+gold^best*da;

% Conditional variance
h=filter([0 a(2)],[1 -a(3)],y*(df-2)/df,s*(df-2)/df)...
  +filter([0 a(1)],[1 -a(3)],ones(t,1));
d=h1.*h2-h.^2;
z=h2.*y1+h1.*y2-2*h.*y;

% Likelihood
if (any(a<low)|any(a>high))
	like=-Inf;
else
  %like=-sum(log(h)+y./h));
  %like=-sum(log(h)+(df+1)*log(1+y./h./df));
	if any(d<=0)|any(1+z./d./df<=0)
		like=-Inf;
	else
		like=-sum(log(d)+(2+df)*log(1+z./d./df))/2;
   end
end

% Gradient
GG=[ filter([0 1],[1 -a(3)],ones(t,1)) ...
     filter([0 1],[1 -a(3)],y*(df-2)/df)	       ...
     filter([0 1],[1 -a(3)],h)	       ];
%g1=y./h.^2-1./h;
%g1=((df+1)*(y./(y+df*h))-1)./h;
 g1=h./d+(2+df)*y./(z+d*df)-(2+df)*h.*z./(z+d*df)./d;
G=GG.*g1(:,ones(1,3));
gra=sum(G);

% Hessian
GG2=GG(:,[1 2 3 1 2 3 1 2 3]).*GG(:,[1 1 1 2 2 2 3 3 3]); 
%g2=-2*y./h.^3+1./h.^2;
%g2=-((df+1)*(y./(y+df*h))-1)./h.^2-(df*(df+1))*(y./(y+df*h).^2./h);
 g2=1./d+2.*h.^2./d.^2-(2+df).*y./(z+d.*df).^2.*(-2.*y-2.*df.*h) ...
   -(2+df).*z./(z+d.*df)./d+2.*(2+df).*h.*y./(z+d.*df)./d ...
   +(2+df).*h.*z./(z+d.*df).^2./d.*(-2.*y-2.*df.*h) ...
   -2.*(2+df).*h.^2.*z./(z+d.*df)./d.^2;
HH=zeros(t,9);
HH(:,3)=filter([0 1],[1 -a(3)],GG(:,1));
HH(:,7)=HH(:,3);
HH(:,6)=filter([0 1],[1 -a(3)],GG(:,2));
HH(:,8)=HH(:,6);
HH(:,9)=filter([0 2],[1 -a(3)],GG(:,3));
H=GG2.*g2(:,ones(1,9))+HH.*g1(:,ones(1,9));
hes=reshape(sum(H),3,3);

% Negative definite
[u,val]=eig(hes);
val=diag(val);
if all(val>0)
   hes=-eye(3);
   negdef=0;
elseif any(val>0)
   negdef=0;
   val=min(val,max(val(val<0))/big);
	hes=u*diag(val)*u';
else
	negdef=1;
end

% Steepest Ascent or Newton
if any([a==low a==high])
   da=-((gra*gra')/(gra*hes*gra'))*gra;
else
   da=-gra/hes;
end

% Termination criterion
term=da*gra';
if ((term<tol1)&negdef)
	break
end

% If you are on the boundary and want to get out, slide along
da((a==low)&(da<0))=zeros(size(da((a==low)&(da<0))));
da((a==high)&(da>0))=zeros(size(da((a==high)&(da>0))));

% If you are stuck in a corner, terminate too
if all(da==0)
   break
end

% Go no further than next boundary
hit=[(low(da~=0)-a(da~=0))./da(da~=0) ...
      (high(da~=0)-a(da~=0))./da(da~=0)];
hit(hit<=0)=[];
da=min([hit 1])*da;

% Step search	
best=0;
newa=a+gold^(best-1)*da;
if (any(newa<low)|any(newa>high))
	left=-Inf;
else
   h=filter([0 newa(2)],[1 -newa(3)],y*(df-2)/df,s*(df-2)/df)...
  	  +filter([0 newa(1)],[1 -newa(3)],ones(t,1));
	d=h1.*h2-h.^2;
	z=h2.*y1+h1.*y2-2*h.*y;
  %left=-sum(log(h)+y./h);
  %left=-sum(log(h)+(df+1)*log(1+y./h./df));
	if any(d<=0)|any(1+z./d./df<=0)
		left=-Inf;
	else
		left=-sum(log(d)+(2+df)*log(1+z./d./df))/2;
   end
end
newa=a+gold^best*da;
if (any(newa<low)|any(newa>high))
	center=-Inf;
else
	h=filter([0 newa(2)],[1 -newa(3)],y*(df-2)/df,s*(df-2)/df)...
	  +filter([0 newa(1)],[1 -newa(3)],ones(t,1));
	d=h1.*h2-h.^2;
	z=h2.*y1+h1.*y2-2*h.*y;
  %center=-sum(log(h)+y./h);
  %center=-sum(log(h)+(df+1)*log(1+y./h./df));
	if any(d<=0)|any(1+z./d./df<=0)
		center=-Inf;
	else
		center=-sum(log(d)+(2+df)*log(1+z./d./df))/2;
   end
end
newa=a+gold^(best+1)*da;
if (any(newa<low)|any(newa>high))
	right=-Inf;
else
	h=filter([0 newa(2)],[1 -newa(3)],y*(df-2)/df,s*(df-2)/df)...
   	  +filter([0 newa(1)],[1 -newa(3)],ones(t,1));
	d=h1.*h2-h.^2;
	z=h2.*y1+h1.*y2-2*h.*y;
  %right=-sum(log(h)+y./h);
  %right=-sum(log(h)+(df+1)*log(1+y./h./df));
	if any(d<=0)|any(1+z./d./df<=0)
		right=-Inf;
	else
		right=-sum(log(d)+(2+df)*log(1+z./d./df))/2;
   end
end
if all(like>[left center right])|all(left>[center right])
	while 1
		best=best-1;
		center=left;
		newa=a+gold^(best-1)*da;
		if (any(newa<low)|any(newa>high))
			left=-Inf;
		else
			h=filter([0 newa(2)],[1 -newa(3)],y*(df-2)/df,s*(df-2)/df)...
		  	  +filter([0 newa(1)],[1 -newa(3)],ones(t,1));
			d=h1.*h2-h.^2;
			z=h2.*y1+h1.*y2-2*h.*y;
        %left=-sum(log(h)+y./h);
        %left=-sum(log(h)+(df+1)*log(1+y./h./df));
			if any(d<=0)|any(1+z./d./df<=0)
				left=-Inf;
			else
				left=-sum(log(d)+(2+df)*log(1+z./d./df))/2;
		   end
		end
		if all(center>=[like left])
			break
		end
	end
elseif all(right>[left center])
	while 1
		best=best+1;
		center=right;
		newa=a+gold^(best+1)*da;
		if (any(newa<low)|any(newa>high))
			right=-Inf;
		else
			h=filter([0 newa(2)],[1 -newa(3)],y*(df-2)/df,s*(df-2)/df)...
		   	  +filter([0 newa(1)],[1 -newa(3)],ones(t,1));
			d=h1.*h2-h.^2;
			z=h2.*y1+h1.*y2-2*h.*y;
        %right=-sum(log(h)+y./h);
        %right=-sum(log(h)+(df+1)*log(1+y./h./df));
			if any(d<=0)|any(1+z./d./df<=0)
				right=-Inf;
			else
				right=-sum(log(d)+(2+df)*log(1+z./d./df))/2;
		   end
		end
		if center>right
			break
		end
	end
end

% If stuck at boundary then stop
%if (center==like)&(any(a<low)|any(a>high))
%	break
%end

% End of optimization loop
%disp(a)
%disp(da)
%disp(gra)
%disp(hes)
%disp([like best iter])
%keyboard

end

q=a;

