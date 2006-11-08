function j = taujacobian(family, alpha)
%
%   FUNCTION J = TAUJACOBIAN(FAMILY, ALPHA)
%   
%   Return the derivative of g'(ALPHA), where
%   where tau = g(ALPHA).
%
%   INPUT
%       FAMILY: One of {'Clayton', 'Frank', 'Gumbel', 'Gaussian', 'AMH',
%              'FGM', 'Arch12', 'Arch14'}.
%       ALPHA:  Copula parameter.
%   
%   OUTPUT
%       J:      The Jacobian evaluated at ALPHA.
%

%   G. Evin, 2005
%   D. Huard, 2006


if nargin < 2
    error('Requires two input arguments');
end

pass = check_alpha(family, alpha);
if any(~pass)
    error('Bad parameters: %s.', mat2str(alpha(~pass)))
end

switch lower(family)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ellipitical copulas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 'gaussian'
        
        j = (2/pi).*(1./sqrt(1-alpha.^2));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Archimedean copulas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'clayton'
        j = 2./((alpha + 2).^2);
        
    case {'frank'} 
        t1 = exp(alpha);
        t2 = (pi^2 + 3*alpha .* ( 1 + alpha + alpha ./ (t1-1)));
        t3 = 6 * alpha .* log(1 - t1);
        t4 = 6 * dilog(t1);
        j =  (4* t2  - t3 - t4)./ 3 ./ alpha.^3;
        
        
    case 'gumbel'
        j = 1./(alpha.^2);
        
    case 'fgm'
        j = 2/9;
        
    case 'amh'
        t1 = alpha.^ 2;
        t3 = log(1-alpha);
        j = -2 / 3 * (t1 + 2 * t3 .* alpha - (2 * alpha) - 2 * t3) ./ t1 ./ alpha;
        
    case 'arch12'
        j = (2/3)./(alpha.^2);
        
    case 'arch14'
        j = 4./((1+2*alpha).^2);   
        
    otherwise
        error('Unrecognized copula type: ''%s''.',family);
end