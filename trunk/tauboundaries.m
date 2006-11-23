function bounds = tauboundaries(family);
%
%   FUNCTION BOUNDS, FUNC = TAUBOUNDARIES(FAMILY)
%
%   Return the minimum and maximum Kendall's tau spanned by a given copula
%   family, and a function handle returning one is tau is in the domain,
%   and 0 otherwise.
% 
%   INPUT
%       FAMILY: A copula from the set {'gaussian', 't', 'clayton', 'gumbel', 'frank', 
%               'arch12', 'arch14', 'joe', 'fgm', 'amh'}
%
%   OUTPUT
%       BOUNDS: [tau_min, tau_max]
%               Domain spanned by Kendall's tau for FAMILY.
%

%   D. Huard, Nov. 2006
%   Pourquoi tau_max de amh etait a .3269 au lieu de 1/3 ?

tau_min = -1;
tau_max = 1;
switch lower(family)    
    case 'amh'
        tau_min = -0.181726;
        tau_max = 1/3;    
    case {'arch12' 'arch14'}
        tau_min = 1/3;
    case 'clayton' 
        tau_min = eps;
    case 'fgm'
        tau_min = -2/9;
        tau_max = 2/9;    
    case 'frank'
        % Discontinuity at 0
        tau_min = -1;
        tau_max = 1;
    case {'gaussian' 't'}        
    case 'gumbel'
        tau_min = 0.;
    otherwise 
        error('Copula family ''%s'' not recognized.', family)
end

bounds =  [tau_min, tau_max];
