function p = posterior(alpha, family, U, prior_tau)
%
%   P = POSTERIOR(FAMILY, U, ALPHA [,PRIOR_TAU])
%
%   Return the posterior density computed at ALPHA.
%   \prod c(u,v|ALPHA) * prior(ALPHA|TAU) * prior(TAU) 
%
%   INPUTS
%       FAMILY   : one of { 'arch12', 'arch14' 'ind', 'fgm' 'gaussian', 'gumbel' 'clayton' 'frank' 'amh'}       
%       U        : Nx2 matrix of quantiles (u,v) in [0,1]^2.
%       ALPHA    : 1xM vector of copula parameters.
%       PRIOR_TAU: Function of TAU returning the normalized prior for TAU.
%                  
%
%   OUTPUT
%       L: Likelihood at each parameter (1xM). 
%

%   G. Evin & D. Huard, 2006

%   Compute density
c = copulapdf(family, U, alpha);
if any(isinf(c))
    error('Inf in copulapdf.')
end
% Compute data likelihood
likelihood = sum(log(c), 1);

% Compute prior for the parameter assuming an uniform prior on TAU.
prior_alpha = log(taujacobian(family, alpha));

% Prior for Tau
pr_tau = log(prior_tau(copulastat(family, alpha)));

% Combine likelihood and priors
p = exp(likelihood + prior_alpha + pr_tau);


