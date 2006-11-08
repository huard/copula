function bounds = tauboundaries(type);
%
% Function bounds = tauboundaries(type)
%
% Return the minimum and maximum Kendall's tau spanned by a given copula family.
% 
% type: A copula from the set:{clayton, gumbel, frank, arch12, arch14, joe, fgm, amh}
%
% Pourquoi tau_max de amh etait a .3269 au lieu de 1/3 ?

tau_min = -1;
tau_max = 1;
switch lower(type)    

    case 'amh'
        tau_min = -0.181726;
        tau_max = 1/3

    case {'arch12' 'arch14'}
        tau_min = 1/3;

    case 'clayton' 
        tau_min = eps;

    case 'fgm'
        tau_min = -2/9;
        tau_max = 2/9;    

    case 'frank'
        tau_min(1) = -1
        tau_min(2) = eps
        tau_max(1) -eps
        tau_max(2) = 1
        
    case 'gaussian'

    case 'gumbel'
	tau_min = 0.

    otherwise 
        error('Copula family not recognized: %s\n', type)

end
