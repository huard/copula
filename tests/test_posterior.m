% TEST POSTERIOR

family = {'gaussian', 'gumbel' 'clayton' 'fgm' 'frank' 'amh' 'arch12', 'arch14'};
a = linspace(-10, 10,200);
N = 50;
for i=1:length(family)
    f = family{i};

    % Select a parameter.
    bounds = tauboundaries(f);
    tau = bounds(1) + .8 * diff(bounds);
    taus = linspace(bounds(1)+.1, bounds(2)-.05, 100);
    alpha = copulaparam(f, taus);
    
    % Generate random values using selected tau.
    U = copularnd(f, copulaparam(f, tau), N);
    
    % Compute likelihood at diffent parameter values.
    like = posterior(alpha, f, U, inline('1'));
    
    % In principle,the likelihood should reach a maximum at the k th value.
    [L, K] = max(like);
    
    if ~isnear(copulastat(f, alpha(K)), tau, .05) 
        fprintf('There may be a problem with copula_like with family ''%s''.\n synth alpha:\t%f\n estimated: \t%f\n', f, alpha(k), alpha(K))
    end
end