function p = bcs(types, data, bounds)
% function probability = bcs(types, data, bounds)
%
% Function that returns the probability of copula models in types, 
% given a vector of bivariate empirical cdf. 
%
% INPUTS:   type is one of the list {'ind' 'gumbel' 'clayton' 'sim' 'frank' 'gb' 'amh' 'joe'}
%           data is the matrix (u,v), the empirical cdf.
%           bounds: integration boundary on Kendall's tau. 
%	     
% OUTPUTS:  The probability of each copula model.


	for i=1:length(types)	
	    [tau_min,tau_max] = int_tau(char(types{i}));
	    parmin = max(bounds(0), tau_min)
	    parmax = min(bounds(1), tau_max)
	    p(i) = quadg('posterior',parmin,parmax,1e-3,[0,128],data,char(types{i}));
	end
