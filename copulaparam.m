function alpha = copulaparam(family,tau)
%
%   ALPHA = COPULAPARAM(FAMILY,TAU)
%   
%   Return the copula parameter, given Kendall's rank correlation.
%
%   INPUT
%       FAMILY: One of {'AMH' 'Arch12' 'Arch14' 'Clayton' 'FGM' 'Frank'
%               'Gaussian' 'Gumbel' 't'}
%       TAU:    Kendall's rank correlation.       
%   
%   OUTPUT
%	    ALPHA:  Copula parameter.       
%   
%
%   Example:
%      % Determine the linear correlation coefficient corresponding to a
%      % bivariate Gaussian copula having a rank correlation of -0.5
%      tau = -0.5
%      rho = copulaparam('gaussian',tau)
%

%   Written by Peter Perkins, The MathWorks, Inc.
%   Revision: 1.0  Date: 2003/09/05
%   This function is not supported by The MathWorks, Inc.
%   Modified by G. Evin & D. Huard, 2006.

if nargin < 2
    error('Requires two input arguments.');
end

pass = check_tau(family, tau);
if any(~pass)
    error('Bad parameters')
end

switch lower(family)
    case {'gaussian' 't'}
        alpha = sin(tau.*pi./2);
        
    case 'clayton'
        alpha = 2*tau ./ (1-tau+eps);
        
    case 'frank'
        for i=1:length(tau)
            if tau(i) == 0
                alpha(i) = 0;
            elseif abs(tau(i)) == 1
                alpha(i) = sign(tau(i)).*realmax;
            else
                % There's no closed form for alpha in terms of tau, so alpha has to be
                % determined numerically.
                alpha(i) = fzero(@invcopulastat,sign(tau(i)),[],'frank', tau(i));
                
            end
        end
        
    case 'fgm'
        alpha = tau.*(9/2);
        
    case 'gumbel'
        alpha = 1 ./ (1-tau);
        
    case {'amh'}
        % There's no closed form for alpha in terms of tau, so alpha has to be
        % determined numerically.
        for i=1:length(tau)
            if tau(i) == 0 
                alpha(i) = 0 ;
            else
                alpha(i) = fzero(@invcopulastat,[-1 0.99999999],[],'amh',tau(i));
            end
        end
               
    case 'arch12'
        alpha = 2./(3.*(1-tau));
        
    case {'arch14'}
        alpha = 1./(1-tau)-1/2;
        
    otherwise
        error('Unrecognized copula family: ''%s''.',family);
end


function err = invcopulastat(alpha, family, target_tau)
% Return difference between target_tau and tau computed from copulastat
% with a guess for alpha.
guess_tau = copulastat(family, alpha);
err = guess_tau - target_tau;


