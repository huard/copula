function c = copulacdf(family, varargin)
%
%   Function C = COPULACDF(FAMILY, U, ALPHA)
%
%   Return  C(u,v|alpha) for the given copula family.
% 
%   INPUT:
%       FAMILY: One from {'AMH' 'Arch12' 'Arch14' 'Clayton' 'FGM' 'Frank'
%               'GB' 'Gumbel' 'Ind' 'Joe'}
%       U: Nx2 vector of quantiles (u,v) in [0,1]^2.
%       ALPHA: 1xM vector of copula parameters
%           or 
%   COPULACDF(FAMILY, U1, U2, ALPHA)
%       U1, U2:  Matrices or row vectors.
%                If U1 is (1xN) and U2 is (1xM), c is (NxM)
%                If U1 is (NxM), U2 must be (N,M), N~=1.
%       ALPHA:   Scalar copula parameter
%
%   OUTPUT
%       C:      Copula CDF C(u,v|alpha) (NxM)
%

%   G. Evin, 2005
%   D. Huard, 2006

% Check alpha is in the domain covered by the family.
if nargin == 3
    U = varargin{1};
    u = U(:,1);
    v = U(:,2);
    alpha = varargin{2};
elseif nargin == 4
    u = varargin{1};
    v = varargin{2};
    alpha = varargin{3};
end

% Check alpha is in the domain covered by the family.
pass = check_alpha(family, alpha);
if ~all(pass)
    error('Some parameters are not valid.\n%f', alpha(~pass))
end

% Check u,v are in [0,1]^2
if any( (u < 0) | (u > 1) | (v < 0) | (v > 1) )
    error('Some quantiles are outside the unit hypercube.')
end

if nargin == 3
    % Shape checking
    [NU, MU] = size(U);
    [NA, MA] = size(alpha);
    
    if MU ~= 2
        error('Bad shape. U is not Nx2, but rather %s.', mat2str(size(U)))
    end
    
    % Reshape ALPHA
    if NA == 1 
        alpha = repmat(alpha, NU, 1);
    elseif NA ~= NU && NU ~= 1
        error('Number of parameters must be 1, identical to number of couples in U, or a row vector.')
    end
    
    % Reshape u,v
    if NU == 1
        u = repmat(u, NA, MA);
        v = repmat(v, NA, MA);
    else
        u = repmat(u, 1, MA);
        v = repmat(v, 1, MA);
    end
elseif nargin == 4
    if ~all(size(alpha)==1)
        error('Alpha must be a scalar.')
    end
    su = size(u);
    sv = size(v);
    if ~any(su==1)
        if all(su ==sv)
            alpha = repmat(alpha, size(u));
        else
            error('If U1 and U2 are matrices, they must have the same size.')
        end
    else
        [u,v] = meshgrid(u,v);
        alpha = repmat(alpha, size(u));
    end
end
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
        % C(u,v) = u^(-alpha)+v^(-alpha)-1)^(-1/alpha)
        c = (u.^(-alpha)+v.^(-alpha)-1).^(-1./alpha);
        
    case 'fgm' % Farlie, Gumbel, Morgenstern
        c = u.*v.*(1+alpha.*(1-u).*(1-v));
        
    case 'amh' % Ali, Mikhail, Haq
        c = (u.*v)./(1-alpha.*(1-u).*(1-v));
        
    case 'gb' % Gumbel & Barnett
        c = u.*v.*exp(-alpha.*log(u).*log(u));
        
    case 'frank' % Frank, 1979
        c = (-1./alpha).*log(1+(exp(-alpha.*u)-1).*(exp(-alpha.*v)-1)./(exp(-alpha)-1));
        
    case 'joe' % Joe, 1997
        c = 1 - (((1-u).^alpha)+((1-v).^alpha)-((1-u).^alpha).*((1-v).^alpha)).^(1./alpha);
        
    case 'arch12' % Archemedean copula #12 in Nelsen's book.
        c = 1 ./ (1 + ((1 ./ u - 1) .^ alpha + (1 ./ v - 1) .^ alpha) .^ (1 ./ alpha));

    case 'arch14' % Archemedean copula #14 in Nelsen's book.
        t1 = 1 ./ alpha;
        t2 = u .^ (-t1);
        t4 = (t2 - 1) .^ alpha;
        t5 = v .^ (-t1);
        t7 = (t5 - 1) .^ alpha;
        t9 = (t4 + t7) .^ t1;
        c = (1 + t9) .^ (-alpha);

    otherwise
        error('Copula family ''%s'' not recognized.', family)
end
