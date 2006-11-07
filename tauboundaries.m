function [tau_min,tau_max] = int_tau(type);
% Function [tau_min, tau_max] = int_tau(type)
%
% Return the minimum and maximum value that can take Kendall's tau for a given copula.
% Input 
%      type: A copula from the set:{clayton, gumbel, frank, arch12, arch14, joe, fgm, amh}

tau_min=-.95;
tau_max=.95;
switch lower(type)
    case {'clayton' 'gumbel'}
        if tau_min < 0
            tau_min = eps;
        end
        
    case 'frank'
        
    case {'arch12' 'arch14'}
        if tau_min < 1/3
            tau_min = 1/3;
        end
        
    case {'joe'}
        if tau_min <= 0
            tau_min = eps;
        end
        
    case 'fgm'
        if tau_min < -2/9
            tau_min = -2/9;
        end
        if tau_max > 2/9
            tau_max = 2/9;
        end
        
    case 'amh'
        if tau_min < -0.1817
            tau_min = -0.1817;
        end
        if tau_max > 0.3269
            tau_max = 0.3269;
        end
end
