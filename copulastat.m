function tau = copulastat(family,alpha)
%
%   TAU = COPULASTAT(FAMILY,ALPHA)
%   
%   Returns Kendall's rank correlation TAU
%
%   INPUT
%       FAMILY: One of 'Clayton', 'Frank', 'Gumbel', 'Gaussian', 't', 'AMH', 'GB', 
%               'Joe', 'FGM', 'Arch12', 'Arch14'.
%       ALPHA: Copula parameter.
%   
%   OUTPUT
%       TAU: Kendall's rank correlation.
%   
%   Example:
%      % Determine the theoretical rank correlation coefficient for a
%      % bivariate Gaussian copula with linear correlation parameter -0.7071
%      rho = -.7071
%      tau = copulastat('gaussian',rho)

%   Original version by Peter Perkins, The MathWorks, Inc.
%   Revision: 1.0  Date: 2003/09/05
%
%   Extended and modified by G.Evin, 2005.
%   Modified by D. Huard, Nov. 2006.
%   Requires MATLAB R13.

if nargin ~= 2
    error('Requires two input arguments.');
end

pass = check_alpha(family, alpha);

if any(~pass)
    error('Invalid parameters.')
end

tau = zeros(size(alpha));

switch lower(family)
    case {'gaussian' 't'}
        tau = 2.*asin(alpha)./pi;
        
    case 'clayton'
        tau = alpha ./ (2 + alpha);
        
    case 'frank'
        i = alpha ~= 0;
        tau(i) = 1 + 4 .* (debye1(alpha(i))-1) ./ alpha(i);
        
    case 'gumbel'
        tau = 1 - 1./alpha;
        
    case 'amh'
for i=1:length(alpha)
        tau(i) = 1 + 4 .* quadg(@lambdaarch,0,1,[],[],'amh', alpha(i));
end
        
    case 'gb'
        for i=1:length(alpha)
	if alpha(i) ~= 0;
	tau = 1 + 4 .* quadg(@lambdaarch,0,1,[],[],'gb', alpha(i));
end
end
	%tau(i) = 1 + 4 .* quad(@(x)lambdaarch('gb', x, alpha(i)),0,1);
        
    case 'joe'
        i = alpha ~= 1;
     for i=1:length(alpha)
	if alpha ~= 1;;
	tau = 1 + 4 .* quadg(@lambdaarch,0,1,[],[],'joe', alpha(i));
end
end

%	tau(i) = 1 + 4 .* quad(@(x)lambdaarch('joe', x, alpha(i)),0,1);

    case 'fgm'    
        tau = (2/9).*alpha;

    case 'arch12'    
        tau = 1-(2/3)./alpha;
        
    case 'arch14'
        tau = 1-2./(1+2.*alpha);

    otherwise
        error('Unrecognized copula family: ''%s''',family);
end