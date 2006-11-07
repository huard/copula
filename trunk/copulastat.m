function tau = copulastat(type,param)
%COPULASTAT Kendall's rank correlation for a copula.
%
%   TAU = COPULASTAT(TYPE,ALPHA) returns Kendall's rank correlation TAU
%   corresponding to a bivariate Archimedean copula with parameter ALPHA.
%   TYPE is one of 'Clayton', 'Frank', 'Gumbel', 'Gaussian', 't', 'AMH', 'GB', 
%   'Joe', 'FGM', 'Arch12', Arch14'.
%
%   Example:
%      % Determine the theoretical rank correlation coefficient for a
%      % bivariate Gaussian copula with linear correlation parameter -0.7071
%      rho = -.7071
%      tau = copulastat('gaussian',rho)
%
%   Written by Peter Perkins, The MathWorks, Inc.
%   Revision: 1.0  Date: 2003/09/05
%   This function is not supported by The MathWorks, Inc.
%
%   Extended and modified by G.Evin, 2005.
%   Requires MATLAB R13.

if nargin < 2
    error('Requires two input arguments.');
end

switch lower(type)
case {'gaussian' 't'}
    rho = param;
    if ((numel(rho) == 1) && (rho < -1 | 1 < rho)) || ((numel(rho) ~= 1) && ~iscor(rho))
        error('RHO must be a correlation coefficient between -1 and 1, or a positive semidefinite correlation matrix.');
    end
    tau = 2.*asin(rho)./pi;
    
case {'clayton' 'frank' 'frank2' 'gumbel' 'amh' 'gb' 'joe' 'fgm' 'arch12' 'arch14'}
    alpha = param;
    switch lower(type)
    case 'clayton'
        if alpha <= 0
            error('ALPHA must be nonnegative for the Clayton copula.');
        end
        tau = alpha ./ (2 + alpha);
    case 'frank'
        if alpha ~= 0
            tau = 1 + 4 .* (debye1(alpha)-1) ./ alpha;
        else
            tau = 0;
        end
    case 'gumbel'
        if alpha < 1
            error('ALPHA must be greater than or equal to 1 for the Gumbel copula.');
        end
        tau = 1 - 1./alpha;
    case 'amh'
        if (alpha < -1)||(alpha>=1)
            error('ALPHA must be in [-1,1[ for the Ali-Mikhail-Haq copula.');
        else
            tau = 1 + 4 .* quadg(@lambdaarch,0,1,[],[],alpha,'amh');
        end
    case 'gb'
        if (alpha < 0)||(alpha>1)
            error('ALPHA must be in [0,1] for the Gumbel-Barnett copula.');
        elseif alpha ~= 0
            tau = 1 + 4 .* quadg(@lambdaarch,0,1,[],[],alpha,'gb');
        else
            tau = 0;
        end
    case 'joe'
        if alpha < 1
            error('ALPHA must be greater than or equal to 1 for the Joe copula.');
        elseif alpha == 1
            tau = 0;
        else 
            tau = 1 + 4 .* quadg(@lambdaarch,0,1,[],[],alpha,'joe');
        end
        
    case 'fgm'
        if abs(alpha) > 1
            error('ALPHA must lie in [-1 1] for the FGM copula.');
        else 
            tau = (2/9).*alpha;
        end  
        
    case 'arch12'
        if alpha < 1
            error('ALPHA must be greater than or equal to 1 for the Arch12 copula.');
        else 
            tau = 1-(2/3)./alpha;
        end
        
    case 'arch14'
        if alpha < 1
            error('ALPHA must be greater than or equal to 1 for the Arch14 copula.');
        else 
            tau = 1-2./(1+2.*alpha);
        end
    end
    
otherwise
    error('Unrecognized copula type: ''%s''',type);
end
