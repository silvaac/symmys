function [Dates,ZeroPrices]=ProcessDB(Yrs2Mat,Start,Stop)
% fixed-income invariance quest from 
% "Risk and Asset Allocation"-Springer (2005), by A. Meucci


load('DB_FixedIncome'); % loads structures of rolling yields to maturity:  Two_Yrs_Rolling, Five_Yrs_Rolling and Ten_Yrs_Rolling   
                        %              with fields:  Dates (vector of double) and Yields (vector of double)

switch Yrs2Mat
    case 2
       RollingStructure=Two_Yrs_Rolling;
    case 5
       RollingStructure=Five_Yrs_Rolling;
    case 10
       RollingStructure=Ten_Yrs_Rolling;
    otherwise
       error('This is impossible');
 end

Times=find(RollingStructure.Dates  >= datenum(Start) & RollingStructure.Dates  <= datenum(Stop));
Yields=RollingStructure.Yields(Times);
Dates=RollingStructure.Dates(Times);

% compute at the given date the price of the bond 
ZeroPrices = exp(-Yields*Yrs2Mat);