function u = copularnd(family, varargin)
%
%    U = COPULARND(FAMILY, ALPHA [, N])
%
%    Generate random quantiles u,v from the given copula family.
% 
%   INPUT
%	    FAMILY: One of {'ind', 'gaussian', 't', 'clayton', 'frank', 'gumbel', 
%		        'fgm', 'amh', 'arch12', 'arch14', 'gb', 'joe'}.
% 	    ALPHA:  Copula parameters (None for independent copula, [rho, nu]
% 	            for the t copula.)
% 	    N:      Number of random variates, defaults to 1.
%
%   OUTPUT
%	    U = [U1, U2]
%
%   EXAMPLES
%	    copularnd('gumbel', 4.5, 10)
%       copularnd('ind', 100)
%       copularnd('t', [.8, 3], 5)
% 

if strcmp(lower(family), 'ind')
    n = varargin{1};
    u = rand(n,2);
else
    % Checking
    alpha = varargin{1};
    pass = check_alpha(family, alpha);
    if any(~pass)
        error('Bad parameter.')
    end
    if numel(alpha) ~= 1 & ~strcmp(lower(family), 't')
        error('ALPHA must be a scalar.');
    end
    
    if nargin >= 3
        n = varargin{2};
    else
        n = 1;
    end
    
    switch lower(family)    
        case 'gaussian'
            u = normcdf(mvnrnd([0 0],[1 alpha; alpha 1],n));
        case 't'
            rho = alpha(1);
            nu = alpha(2);
            u = tcdf(mvtrnd([1 rho; rho 1],nu,n),nu);
        otherwise
            % Archimedean copulas
            %
            % Random pairs from these copulae can be generated sequentially: first
            % generate u1 as a uniform r.v.  Then generate u2 from the conditional
            % distribution F(u2 | u1; alpha) by generating uniform random values, then
            % inverting the conditional CDF.
            u1 = rand(n,1);
            p = rand(n,1);
            switch lower(family)
                case 'clayton'
                    % The inverse conditional CDF has a closed form for this
                    % copula.
                    if alpha > sqrt(eps)
                        u2 = u1.*(p.^(-alpha./(1+alpha)) - 1 + u1.^alpha).^(-1./alpha);
                    else
                        u2 = p;
                    end
                    u = [u1 u2];
                    
                case 'frank'
                    % The inverse conditional CDF has a closed form for this
                    % copula.
                    if abs(alpha) > log(realmax)
                        u2 = (u1 < 0) + sign(alpha).*u1; % u1 or 1-u1
                    elseif abs(alpha) > sqrt(eps)
                        u2 = -log((exp(-alpha.*u1).*(1-p)./p + exp(-alpha))./(1 + exp(-alpha.*u1).*(1-p)./p))./alpha;
                    else
                        u2 = p;
                    end
                    u = [u1 u2];
                    
                case {'gumbel' 'fgm' 'amh' 'arch12' 'arch14' 'gb' 'joe'}
                    % The inverse conditional CDF does not have a closed form for this
                    % copula.  The inversion must be done numerically.
                    u2 = condCDFinv(@conditionalcdf,u1,p,alpha,family);
                    u = [u1 u2]; 
            end           
    end
end


function u2 = condCDFinv(condCDF,u1,p,par, family)
%
%    CONDCDFINV Inverse conditional distribution function
%    
%    U2 = CONDCDFINV(CONDCDF,U1,P,ALPHA) returns U2 such that
%
%      CONDCDF(U1,U2,ALPHA) = P,
%
%  where CONDCDF is a function handle to a function that computes the
%  conditional cumulative distribution function of U2 given U1, for an
%  archimedean copula with parameter ALPHA.
%
% CONDCDFINV uses a simple binary chop search.  Newton's method or the
% secant method would probably be faster.

lower = zeros(size(p));
upper = ones(size(p));
width = 1;
tol = 1e-12;
while width > tol
    u2 = .5 .* (lower + upper);
    lo = feval(condCDF, family, u1,u2,par) < p;
    lower(lo) = u2(lo);
    upper(~lo) = u2(~lo);
    width = .5 .* width;
end
