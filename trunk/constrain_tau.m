function  newbounds = constrain_tau(family, bounds)
%
%   FUNCTION NEWBOUNDS = CONSTRAIN_TAU(FAMILY, BOUNDS)
%
%   Return boundaries on tau inside the domain range of each copula family.
%
%   INPUT
%       FAMILY: One of {'AMH' 'Arch12' 'Arch14' 'Clayton' 'FGM' 'Frank'
%               'Gaussian' 'Gumbel' 't'}
%       BOUNDS: [tau1, tau2]
%
%   OUTPUT
%       NEWBOUNDS: Intersection of BOUNDS and the FAMILY domain in tau.
%
%   Example:
%       constrain_tau('fgm', [0,1])
%       ans = 
%           0.0000, .22222
%

% D. Huard, Nov. 2006

boundaries = tauboundaries(family);
tau1 = boundaries(1);
tau2 = boundaries(2);
tau_min = max(bounds(1), tau1);
tau_max = min(bounds(2), tau2);
newbounds = [tau_min, tau_max];
