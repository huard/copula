function p = copula_like(family, U, alpha)
%
%   P = COPULA_LIKE(FAMILY, U, ALPHA)
%
%   Return the likelihoods computed at ALPHA.
%   \prod c(u,v|ALPHA) * prior(ALPHA) 
%
%   INPUTS
%       FAMILY: one of {'ind', 'gaussian', 'gumbel' 'clayton' 'sim' 'frank' 'gb' 'amh' 'joe'}       
%       U: Nx2 vector (u,v) in [0,1]^2.
%       ALPHA: 1xM vector of copula alphaameter.
%
%   OUTPUTS
%       P: Vector (1xM) 
%

%   TODO: Implement user-defined prior on tau.
%   G. Evin & D. Huard, 2006

%   Compute density
c = copulapdf(family, U, alpha)

% Compute data likelihood
likelihood = sum(log(c), 1)

% Compute prior on the parameter assuming an uniform prior on TAU.
prior_alpha = log(taujacobian(family, alpha))

% Prior on Tau
% Implement user-defined prior on tau.
prior_tau = log(ones(alpha))

% Combine
loglike = likelihood + prior_alpha + prior_tau

p = exp(loglike);
