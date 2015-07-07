function e=FitError(gamma,X,Y,Correlation);

e = sum(trace(  (exp(-gamma*abs(X-Y))-Correlation)*(exp(-gamma*abs(X-Y))-Correlation)'  ));

