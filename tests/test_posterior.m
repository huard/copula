% TEST POSTERIOR
% Generate random values with some synthetic parameter and check that the
% parameter estimated from the sample correponds well to the synthetic parameter. 

family = {'gaussian', 'gumbel' 'clayton' 'fgm' 'frank' 'amh' 'arch12', 'arch14'};
a = linspace(-10, 10,200);
N = 200;
for i=1:length(family)
    f = family{i};

    % Select a synthetic parameter.
    bounds = tauboundaries(f);
    tau = bounds(1) + .8 * diff(bounds);
    alpha = copulaparam(f, tau);
    
    % Generate random values using the synthetic parameter.
    U = copularnd(f, alpha, N);
    
    % Compute likelihood at diffent parameter values.
    taus = linspace(bounds(1)+.1, bounds(2)-.05, 100);
    alphas = copulaparam(f, taus);
    like = posterior(alphas, f, U, inline('1'));
    
    % In principle,the maximum likelihood should be reached near the
    % synthetic parameter.
    [L, K] = max(like);
    
    if ~isnear(copulastat(f, alphas(K)), tau, .05) 
        fprintf('There may be a problem with posterior.m with family ''%s''.\n synth alpha:\t%f (%.2f)\n estimated: \t%f (%.2f)\n', f, alpha, tau, alphas(K), taus(K))
    end
end