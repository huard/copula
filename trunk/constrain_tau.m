function  newbounds = constrain_tau(family, bounds)
%
%   FUNCTION NEWBOUNDS = CONSTRAIN_TAU(FAMILY, BOUNDS)
%
%   Return boundaries on tau inside the domain range of each copula family.
%
%   INPUT
%       FAMILY: One of the set 
%       BOUNDS: [tau1, tau2]
%
%   OUTPUT
%       NEWBOUNDS: Intersection of BOUNDS and the FAMILY domain in tau.
%
%   Example:
%       restrict_tau('fgm', [0,1])
%       ans = 
%           0.0000, .22222
%

% D. Huard, Nov. 2006

	[tau1, tau2] = tauboundaries(family);
	tau_min = max(bounds(1), tau1);
	tau_max = min(bounds(2), tau2);
	newbounds = [tau_min, tau_max];
