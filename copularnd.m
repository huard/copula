function u = copularnd(family,varargin)
%
%    U = COPULARND(FAMILY, ALPHA, N)
%
%    Generate random u,v from copula family.
% 
%    INPUT
%	FAMILY: One of {'ind', 'gaussian', 't', 'clayton', 'frank', 'gumbel', 
%		'fgm', 'amh', 'arch12', 'arch14', 'gb', 'joe'}.
% 	ALPHA:  Copula parameters (None for independent copula).
% 	N:      Number of random variates.
%
%    OUTPUT
%	U = [U1, U2]
%
%    EXAMPLE
%	copularnd('gumbel', 4.5, 10)
% 

switch lower(family)
    case 'ind'
        n = varargin{1};
        u = rand(n,2);

    case 'gaussian'
        if nargin < 3
            error('Requires three input arguments for the Gaussian copula.');
        end
        rho = varargin{1};
        n = varargin{2};
        if numel(rho) == 1
            if (rho < -1 | 1 < rho)
                error('RHO must be a correlation coefficient between -1 and 1, or a positive semidefinite correlation matrix.');
            end
            u = normcdf(mvnrnd([0 0],[1 rho; rho 1],n));
        else
            if ~iscor(rho)
                error('RHO must be a correlation coefficient between -1 and 1, or a positive semidefinite correlation matrix.');
            end
            p = size(rho,1);
            u = normcdf(mvnrnd(zeros(1,p),rho,n));
        end

    case 't'
        if nargin < 4
            error('Requires four input arguments for the t copula.');
        end
        rho = varargin{1};
        nu = varargin{2};
        n = varargin{3};
        if ~(0 < nu)
            error('NU must be positive for the t copula.');
        end
        if numel(rho) == 1
            if (rho < -1 | 1 < rho)
                error('RHO must be a correlation coefficient between -1 and 1, or a positive semidefinite correlation matrix.');
            end
            u = tcdf(mvtrnd([1 rho; rho 1],nu,n),nu);
        else
            if ~iscor(rho)
                error('RHO must be a correlation coefficient between -1 and 1, or a positive semidefinite correlation matrix.');
            end
            u = tcdf(mvtrnd(rho,nu,n),nu);
        end

        % Archimedean copulas
        %
        % Random pairs from these copulae can be generated sequentially: first
        % generate u1 as a uniform r.v.  Then generate u2 from the conditional
        % distribution F(u2 | u1; alpha) by generating uniform random values, then
        % inverting the conditional CDF.
    case {'clayton' 'frank' 'gumbel' 'fgm' 'amh' 'arch12' 'arch14' 'gb' 'joe'}
        if nargin < 3
            error('Requires three input arguments for an Archimedean copula.');
        end
        alpha = varargin{1};
        if numel(alpha) ~= 1
            error('ALPHA must be a scalar.');
        end

        try indicatrice(copulastat(family,alpha),family);
        catch err = lasterror;
            display(err.message);
            display('The copula parameter is wrong in function copularnd ');
            return;
        end
        n = varargin{2};
        u1 = rand(n,1);

        switch lower(family)
            case 'clayton'
                % The inverse conditional CDF has a closed form for this copula.
                p = rand(n,1);
                if alpha > sqrt(eps)
                    u2 = u1.*(p.^(-alpha./(1+alpha)) - 1 + u1.^alpha).^(-1./alpha);
                else
                    u2 = p;
                end
                u = [u1 u2];

            case 'frank'
                % The inverse conditional CDF has a closed form for this copula.
                p = rand(n,1);
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
                p = rand(n,1);
                u2 = condCDFinv(@copCondCDF,u1,p,alpha,family);
                u = [u1 u2];
        end

    otherwise
        error('Unrecognized copula family: ''%s''',family);
end
