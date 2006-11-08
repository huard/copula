function param = copulaparam(type,tau)
%COPULAPARAM Copula parameter, given Kendall's rank correlation.
%   RHO = COPULAPARAM('Gaussian',TAU) returns the linear correlation
%   parameter RHO corresponding to a Gaussian copula having Kendall's rank
%   correlation TAU.  If TAU is a scalar correlation coefficient, RHO is a
%   scalar correlation coefficient corresponding to a bivariate copula.  If
%   TAU is a P-by-P correlation matrix, RHO is a P-by-P correlation matrix
%   corresponding to a P-variate copula.
%
%   RHO = COPULAPARAM('t',TAU) returns the linear correlation parameter RHO
%   corresponding to a t copula having Kendall's rank correlation TAU.  If
%   TAU is a scalar correlation coefficient, RHO is a scalar correlation
%   coefficient corresponding to a bivariate copula.  If TAU is a P-by-P
%   correlation matrix, RHO is a P-by-P correlation matrix corresponding to
%   a P-variate copula.
%   
%   ALPHA = COPULAPARAM(TYPE,TAU) returns the copula parameter ALPHA
%   corresponding to a bivariate Archimedean copula having Kendall's rank
%   correlation TAU.  TYPE is one of 'Clayton', 'Frank', or 'Gumbel'.
%
%   Example:
%      % Determine the linear correlation coefficient corresponding to a
%      % bivariate Gaussian copula having a rank correlation of -0.5
%      tau = -0.5
%      rho = copulaparam('gaussian',tau)
%
%      % Generate dependent beta random values using that copula
%      u = copularnd('gaussian',rho,100)
%      b = betainv(u,2,2)
%
%      % Verify that those pairs have a sample rank correlation approximately
%      % equal to tau
%      tau_sample = kendall(b)

%   Written by Peter Perkins, The MathWorks, Inc.
%   Revision: 1.0  Date: 2003/09/05
%   This function is not supported by The MathWorks, Inc.
%
%   Requires MATLAB R13.

if nargin < 2
    error('Requires two input arguments.');
end

switch lower(type)
    case {'gaussian' 't'}
        if ((numel(tau) == 1) && (tau < -1 | 1 < tau)) || ((numel(tau) ~= 1) && ~iscor(tau))
            error('TAU must be a correlation coefficient between -1 and 1, or a positive semidefinite correlation matrix.');
        end
        param = sin(tau.*pi./2);
        
    case {'clayton' 'frank' 'fgm' 'gumbel' 'amh' 'gb' 'joe'  'arch12' 'arch14_2' 'arch14'}
        if (tau < -1 | 1 < tau)
            error('TAU must be a correlation coefficient between -1 and 1.');
        end
        switch lower(type)
            case 'clayton'
                nn = (tau == 0) ;
                param(nn) = 0;
                param(~nn) = 2*tau(~nn) ./ (1-tau(~nn));
            case 'frank'
                nn = (tau == 0) ;
                param(nn) = 0;
                if abs(tau) < 1
                    % There's no closed form for alpha in terms of tau, so alpha has to be
                    % determined numerically.
                    warn = warning('off','MATLAB:fzero:UndeterminedSyntax');
                    param = fzero(@frankRootFun,sign(tau(~nn)),[],tau);
                    warning(warn);
                else
                    param(~nn) = sign(tau(~nn)).*Inf;
                end
            case 'fgm'
                if (tau < -2/9) || (tau > 2/9)
                    warning('TAU must lie in [-2/9,2/9] for the FGM copula.');
                end
                param = tau.*(9/2);    
            case 'gumbel'
                %                 if tau < 0
                %                     warning('TAU must be nonnegative for the Gumbel copula.');
                %                     param = 0;
                %                 end
                param = 1 ./ (1-tau);
            case {'amh'}
                if (tau < -0.1817) || (tau > 1/3)
                    error('TAU must lie in [-0.1817,1/3] for the Ali-Mikhail-Haq copula.'); 
                end
                
                nn = (tau == 0) ;
                param(nn) = 0;
                
                % There's no closed form for alpha in terms of tau, so alpha has to be
                % determined numerically.
                warn = warning('off','MATLAB:fzero:UndeterminedSyntax');
                param(~nn) = fzero(@archRootFun,[-1 0.99999999],[],'amh',tau(~nn));
                warning(warn);
            case {'gb'}
                if (tau < -0.3613) || (tau > 0)
                    warning('TAU must lie in [-0.3613,0] for the Gumbel-Barnett copula.'); 
                end
                
                nn = (tau == 0) ;
                param(nn) = 0;
                
                % There's no closed form for alpha in terms of tau, so alpha has to be
                % determined numerically.
                warn = warning('off','MATLAB:fzero:UndeterminedSyntax');
                param(~nn) = fzero(@archRootFun,[0 1],[],'gb',tau(~nn));
                warning(warn);
            case {'joe'}
                if tau < 0 
                    warning('TAU must be nonnegative for the Joe copula.'); 
                end
                
                nn = (tau == 0) ;
                param(nn) = 0;
                % There's no closed form for alpha in terms of tau, so alpha has to be
                % determined numerically.
                warn = warning('off','MATLAB:fzero:UndeterminedSyntax');
                param(~nn) = fzero(@archRootFun,[1 99],[],'joe',tau(~nn));
                % The interval of the parameter is limited to 99 , because of the computation
                % time, alpha equals to 99 corresponds to a tau equals
                % to  0.9895
                warning(warn);
            case 'arch12'
                if tau < 1/3
                    error('TAU must be greater than or equal to 1/3 for the arch12 copula.');
                else
                    param = 2./(3.*(1-tau));
                end
                
                
            case {'arch14_2' 'arch14'} 
                if tau < 1/3
                    error('TAU must be greater than or equal to 1/3 for the arch14 copula.');
                else
                    param = 1./(1-tau)-1/2;
                end 
                
            otherwise
                error('Unrecognized copula type: ''%s''',type);
        end
end