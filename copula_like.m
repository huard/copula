function p = copula_like(alpha, family, U, prior_tau_f)
%
%   P = COPULA_LIKE(FAMILY, U, ALPHA)
%
%   Return the likelihoods computed at ALPHA.
%   \prod c(u,v|ALPHA) * prior(ALPHA|TAU) * prior(TAU) 
%
%   INPUTS
%       FAMILY   : one of {'ind', 'gaussian', 'gumbel' 'clayton' 'sim' 'frank' 'gb' 'amh' 'joe'}       
%       U        : Nx2 matrix of quantiles (u,v) in [0,1]^2.
%       ALPHA    : 1xM vector of copula parameters.
%       PRIOR_TAU: Function of TAU returning the normalized prior for TAU.
%
%   OUTPUT
%       P: Vector (1xM) 
%

%   TODO: Implement user-defined prior on tau.
%   G. Evin & D. Huard, 2006

%   Compute density
c = copulapdf(family, U, alpha);

% Compute data likelihood
likelihood = sum(log(c), 1);

% Compute prior for the parameter assuming an uniform prior on TAU.
prior_alpha = log(taujacobian(family, alpha));

% Prior for Tau
if nargin < 4
    prior
prior_tau = log(prior_tau(copulastat(family, alpha)))

% Combine likelihood and priors
loglike = likelihood + prior_alpha + prior_tau;

p = exp(loglike);
