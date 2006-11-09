function c = copulacdf(family, U, alpha)
%
%   Function C = COPULACDF(FAMILY, U, ALPHA)
%
%   Return  C(u,v|alpha) for the given copula family.
% 
%   INPUT:
%       FAMILY: One of 'Ind', 'Gumbel', 'Clayton', 'FGM', 'AMH', 'GB', 'Frank', 'Joe'.
%       U:      Quantiles (Nx2)       
%       ALPHA:  Copula parameter (1xM)
% 
%   OUTPUT
%       C:      Copula CDF (NxM)
%

%   G. Evin, 2005
%   D. Huard, 2006

% Check alpha is in the domain covered by the family.
pass = check_alpha(family, alpha);
if ~all(pass)
    error('Some parameters are not included in the domain.\n%f', alpha(~pass))
end

% Check u,v are in [0,1]^2
if any( (U < 0) | (U > 1) )
    error('Some quantiles are outside the unit hypercube.')
end

sU = size(U);
sA = size(alpha);

N = sU(1);
M = sA(2);
L = sA(1);

if L > 1 && L ~= N
   error('Number of parameters must be 1, identical to number of couples in U, or a row vector.')
end

% Shape vectors into matrices to compute the pdf for all values of (u,v)
% and all parameters without loops.
u = U(:,1);
v = U(:,2);
u = repmat(u, 1, M);
v = repmat(v, 1, M);
alpha = repmat(alpha, N, 1);

% Compute copula cdf
switch lower(family)
    case 'ind' % Independent
        c = u.*v;
        
    case 'gumbel' % Gumbel       
        a1 = (-log(u)).^alpha;
        a2 = (-log(v)).^alpha;
        a3 = 1./alpha;
        c = exp(-(a1+a2).^a3);
        
    case 'clayton'
        c = (u.^(-alpha)+v.^(-alpha)-1).^(-1./alpha);
        
    case 'fgm' % Farlie, Gumbel, Morgenstern
        c = u.*v.*(1+alpha.*(1-u).*(v));
        
    case 'amh' % Ali, Mikhail, Haq
        c = (u.*v)./(1-alpha.*(1-u).*(v));
        
    case 'gb' % Gumbel & Barnett
        c = u.*v.*exp(-alpha.*log(u).*log(u));
        
    case 'frank' % Frank, 1979
        c = (-1./alpha).*log(1+(exp(-alpha.*u)-1).*(exp(-alpha.*v)-1)./(exp(-alpha)-1));
        
    case 'joe' % Joe, 1997
        c = 1 - (((1-u).^alpha)+((1-v).^alpha)-((1-u).^alpha).*((1-v).^alpha)).^(1./alpha);
        
    otherwise
        error('Copula family ''%s'' not recognized.', family)
end
