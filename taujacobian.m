function j = taujacobian(family, alpha)
%
%   FUNCTION J = TAUJACOBIAN(FAMILY, ALPHA)
%   
%   Return the derivative of tau(alpha) with respect to alpha.
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
%         t1 = exp(alpha);
%         t2 = (pi^2 + 3*alpha .* ( 1 + alpha + alpha ./ (t1-1)));
%         t3 = 6 * alpha .* log(1 - t1);
%         t4 = 6 * dilog(t1);
%         j =  (4* t2  - t3 - t4)./ 3 ./ alpha.^3;
        
t2 = exp(alpha);
t5 = pi ^ 2;
t7 = 1 - t2;
t8 = log(t7);
t9 = alpha .* t8;
t13 = dilog(t2);
t17 = alpha .^ 2;
j = -4 / 3 * (-3 * alpha + 3 * alpha .* t2 - t5 + t5 .* t2 + 6 * t9 - 6 * t9 .* t2 + 6 * t13 - 6 * t13 .* t2 + 3 * t2 .* t17) ./ t17 ./ alpha ./ t7;

        
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

if any(j<0)
    error('Negative value in taujacobian')
end