function alpha = copulaparam(family,tau)
%
%   ALPHA = COPULAPARAM(FAMILY,TAU)
%   
%   Returns copula parameter, given Kendall's rank correlation.
%
%   INPUT
%       FAMILY: One of 'Clayton', 'Frank', 'Gumbel', 'Gaussian', 't', 'AMH', 'GB', 
%               'Joe', 'FGM', 'Arch12', 'Arch14'.
%       TAU: Kendall's rank correlation.%       
%   
%   OUTPUT
%	ALPHA: Copula parameter.       
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

pass = check_tau(family, tau)
if any(~pass)
    error('Bad parameters')
end

switch lower(family)
    case {'gaussian' 't'}
        alpha = sin(tau.*pi./2);
  
    case 'clayton'
        alpha = 2*tau ./ (1-tau);
                
            case 'frank'
                nn = (tau == 0) ;
                alpha(nn) = 0;
                if abs(tau) < 1
                    % There's no closed form for alpha in terms of tau, so alpha has to be
                    % determined numerically.
                    warn = warning('off','MATLAB:fzero:UndeterminedSyntax');
                    alpha = fzero(@frankRootFun,sign(tau(~nn)),[],tau);
                    warning(warn);
                else
                    alpha(~nn) = sign(tau(~nn)).*Inf;
                end
                
            case 'fgm'
                alpha = tau.*(9/2);
                
            case 'gumbel'
                alpha = 1 ./ (1-tau);
                
            case {'amh'}
                nn = (tau == 0) ;
                alpha(nn) = 0;
                % There's no closed form for alpha in terms of tau, so alpha has to be
                % determined numerically.
                warn = warning('off','MATLAB:fzero:UndeterminedSyntax');
                alpha(~nn) = fzero(@archRootFun,[-1 0.99999999],[],'amh',tau(~nn));
                warning(warn);
                
            case {'gb'}
                nn = (tau == 0) ;
                alpha(nn) = 0;
                % There's no closed form for alpha in terms of tau, so alpha has to be
                % determined numerically.
                warn = warning('off','MATLAB:fzero:UndeterminedSyntax');
                alpha(~nn) = fzero(@archRootFun,[0 1],[],'gb',tau(~nn));
                warning(warn);
                
            case {'joe'}
                nn = (tau == 0) ;
                alpha(nn) = 0;
                % There's no closed form for alpha in terms of tau, so alpha has to be
                % determined numerically.
                warn = warning('off','MATLAB:fzero:UndeterminedSyntax');
                alpha(~nn) = fzero(@archRootFun,[1 99],[],'joe',tau(~nn));
                % limited to 99 for computation limit. alpha = 99
                % corresponds to tau = 0.9895
                warning(warn);
                
            case 'arch12'
                    alpha = 2./(3.*(1-tau));
 
            case {'arch14'}
                    alpha = 1./(1-tau)-1/2;

            otherwise
                error('Unrecognized copula family: ''%s''',family);
        end
end
