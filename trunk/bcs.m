function p = bcs(family, U, boundaries, prior_tau)
%
%   FUNCTION BCS(FAMILY, U, BOUNDARIES, PRIOR_TAU)
%
%   Bayesian Copula Selection
%   Return the weight of each copula family.
%
%   INPUT   
%       FAMILY    : Subset of {'Ind' 'Gumbel' 'Clayton' 'Frank' 'AMH' 'Joe' 'Arch12' 'Arch14'}
%       U         : Nx2 matrix of quantiles (u,v).
%       BOUNDARIES: Integration boundaries on Kendall's tau 
%                   Default: [-.95,.95])
%       PRIOR_TAU : Function handle returning the parent prior for TAU. 
%                   Default is uniform over the domain of tau covered by each copula.
%
%   OUTPUT
%       P         : The weight of each copula family.
%
%   NOTES
%       A prior for the model is applied, proportional to 1/Lambda,
%       where Lambda is the interval spanned by Kendall's tau. This 
%       interval is equal the intersection of the copula boundaries 
%       on tau and the boundaries specified by the user. 
%
%       For some data and parameters, the pdf returned by copulapdf is Inf
%       or 0, which ruins the integration. To avoid this problem, try to
%       restrict the integration domain (do not use [-1,1], this won't
%       work, [.-97, .97] is a better choice, and even then, problems may
%       appear).

%   Reference
%   Huard, D., Évin, G. and Favre, A-C. Bayesian Copula Selection, 
%   Journal of Computational Statistics and Data Analysis, 2005, 51, 809-822.


% Set defaults
if nargin <= 3
    prior_tau = inline('1');
end

if nargin <= 2
    boundaries = [-.95, .95];
end

% Make sure U is in the unit hypercube.
if any( (U < 0) | (U > 1) )
    error('Some quantiles are outside the unit hypercube.')
end

% Make sure we can iterate on family
if ~iscell(family)
    family = {family}
end

% Loop on each family to compute individual weight.
for i=1:length(family)
    if strcmp(lower(family(i)), 'ind')
        p(i) = 1;
    else
        % Constrain the boundaries to the domain covered by each family.
        bounds_tau = constrain_tau(family{i}, boundaries);
        
        % If bounds_tau is contrained, shift the boundaries away from the limit
        % of the domain, since many functions don't deal well with those.
        shift = bounds_tau ~= boundaries;
        correct = [1,-1];
        bounds_tau(shift) = bounds_tau(shift) .* (1 + correct(shift).*sign(bounds_tau(shift))*.01);
        
        % Translate the boundaries on tau in copula parameters.
        bounds_alpha = copulaparam(family{i}, bounds_tau);
        alpha_min = bounds_alpha(1);
        alpha_max = bounds_alpha(2);
        
        % Integrate the likelihood over the parameter range.
        p(i) = quadg('posterior',alpha_min, alpha_max, 1e-4, [0,128], family{i}, U, prior_tau);
        
        % Prior for the copula family
        p(i) = p(i)/diff(bounds_tau);
    end
end
