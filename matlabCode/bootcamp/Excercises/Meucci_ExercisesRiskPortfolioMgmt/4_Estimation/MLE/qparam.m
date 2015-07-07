function q=qparam(p,th)

if th<=0
    m=th;
    s=sqrt(th^2);
    nu=1;
    q=m+s*tinv(p,nu);
else
    m=th;
    s=sqrt((th-.01)^2);
    q=logninv(p,m,s);
end