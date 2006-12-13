%TEST_BCS

% Check that the method identifies correctly the copula that generated the
% random sample.

family_generated = {'gaussian', 'gumbel' 'clayton' 'fgm' 'frank' 'amh' 'arch12', 'arch14'};
family_tested = {'gaussian', 'gumbel' 'clayton' 'fgm' 'frank' 'amh' 'arch12', 'arch14', 'ind'};

% Number of trials
T=5;

% Number of samples
N = 100;

warning off MATLAB:fzero:UndeterminedSyntax
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
        p(i, :, k)  = bcs(family_tested, U, [0.01, .97]);
        P(i,:,k) = p(i,:,k)/sum(p(i,:,k));
    end
end

success = zeros(length(family_generated),1);
% Count the ratio of successful identification
for i=1:length(family_generated)
    for k=1:T
        [Y,I] = max(p(i, :, k));
        if (I == i)
            success(i)  = success(i) + 1;
        end
    end
end

fprintf('%-19s%10s\n', 'Family', 'Number of successful identifications')
fprintf('-------------------------------------------------------\n')
for i=1:length(family_generated)
    fprintf('%-10s%10d/%d\n', family_generated{i}, success(i), T)
end