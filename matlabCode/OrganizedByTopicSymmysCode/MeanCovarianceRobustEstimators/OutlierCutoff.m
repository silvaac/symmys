function Omega = OutlierCutoff(d,d_0)

IndexFlat=find(d<=d_0);

b_2=1.25;
Omega = (d_0./d) .* exp(-.5*((d-d_0).^2)/(b_2^2)  );

Omega(IndexFlat)=1;


