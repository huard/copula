function  newbounds = restrict_tau(type, bounds)
%
% newbounds = restrict_tau(bounds)
%
% Return boundaries on tau inside the domain range of each copula family.
%
% Example:
% restrict_tau('fgm', [0,1])
% ans = 
%   0.0000, .22222
%

% D. Huard, Nov. 2006

	[tau1, tau2] = tauboundaries(type);
	tau_min = max(bounds(1), tau1);
	tau_max = min(bounds(2), tau2);
	newbounds = [tau_min, tau_max];
