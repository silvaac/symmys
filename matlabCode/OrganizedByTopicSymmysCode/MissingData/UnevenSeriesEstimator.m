function [E_hat,V_hat]=UnevenSeriesEstimator(Z);
% this function estimates the parameters of a multivariate normal distribution 
% from an unbalanced panel of time series of different length
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% re-order panel
[T,N]=size(Z);
Num_Nan=[];
for n=1:N
    Num_Nan=[Num_Nan sum(0+isnan(Z(:,n)))];
end
[Miss,Index]=sort(Num_Nan);

R=Z(:,Index);
j=1;
Panel(1).Index =1;
Panel(1).s=1;
for n=2:N
    if Miss(n)==Miss(n-1)
        Panel(j).Index=[Panel(j).Index n];
    else
        j=j+1;
        Panel(j).Index=n;
        Panel(j).s=Miss(n)+1;
    end
end
J=j;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate
Panel(1).E_hat_sqb=mean(R(Panel(1).s:end,Panel(1).Index))';
Panel(1).V_hat_sqb=cov(R(Panel(1).s:end,Panel(1).Index),1);

for j=2:J

    Y = R(Panel(j).s:end,Panel(j).Index);
    X = [ones(T-Panel(j).s+1,1) R(Panel(j).s:end,1:Panel(j-1).Index(end))];
    C_hat=inv(X'*X)*X'*Y;
    a_hat=C_hat(1,:)';
    B_hat=C_hat(2:end,:)';
    Sigma_hat=1/(T-Panel(j).s+1)*(Y-X*C_hat)'*(Y-X*C_hat);
    E_hat=a_hat+B_hat*Panel(j-1).E_hat_sqb;

    Panel(j).E_hat_sqb=[Panel(j-1).E_hat_sqb
        E_hat];
    Panel(j).V_hat_sqb=[Panel(j-1).V_hat_sqb  Panel(j-1).V_hat_sqb*B_hat'
        B_hat*Panel(j-1).V_hat_sqb    Sigma_hat+B_hat*Panel(j-1).V_hat_sqb*B_hat'];

end

E_hat=Panel(J).E_hat_sqb;
V_hat=Panel(J).V_hat_sqb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% re-order panel back
for n=1:N
    CounterIndex(n)=find(Index==n);
end

E_hat=E_hat(CounterIndex);
V_hat=V_hat(CounterIndex,CounterIndex);