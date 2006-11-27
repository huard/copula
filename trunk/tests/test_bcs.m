%TEST_BCS

% Check that the method identifies correctly the copula that generated the
% random sample.

family_generated = {'gaussian', 'gumbel' 'clayton' 'fgm' 'frank' 'amh' 'arch12', 'arch14'};
family_tested = {'gaussian', 'gumbel' 'clayton' 'fgm' 'frank' 'amh' 'arch12', 'arch14', 'ind'};

% Number of trials
T=1;

% Number of samples
N = 100;

for i=1:length(family_generated)
    f = family_generated{i};
    
    % Select a synthetic parameter.
    bounds = tauboundaries(f);
    tau = bounds(1) + .8 * diff(bounds);
    alpha = copulaparam(f, tau);
    
    for k=1:T
        % Generate random values using the synthetic parameter.
        U = copularnd(f, alpha, N);
        
        % Compute weight of each copula family.
        p(i, :, k)  = bcs(family_tested, U, [0.01, .97])
    end
end

