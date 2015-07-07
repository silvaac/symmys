function Allocation = DecisionBestPerfomer(Market,InvestorProfile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find index of best performer
[a,B]=max(Market.LinRetsSeries(end,:));

% invest in that asset
I=eye(length(Market.CurrentPrices));
Allocation = InvestorProfile.Budget*I(:,B)/Market.CurrentPrices(B);





