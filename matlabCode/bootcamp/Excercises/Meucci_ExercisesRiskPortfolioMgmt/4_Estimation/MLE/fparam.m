function f=fparam(x,th)

if th<=0
    m=th;
    s=sqrt(th^2);
    nu=1;
    f=1/s*tpdf((x-m)/s,nu);
else
    m=th;
    s=sqrt((th-.01)^2);
    f=lognpdf(x,m,s);
end