% Test copularnd

family = {'gaussian', 'clayton', 'frank', 'gumbel',... 
		        'fgm', 'amh', 'arch12', 'arch14'};

%Test correlation of generated sample.
%It should match the tau corresponding to the parameter given to
%copularnd.
N = 800;
for i=1:length(family)
    f = family{i};
    bounds= tauboundaries(f);
    tau = bounds(1) + .9 * diff(bounds);
    alpha = copulaparam(f, tau);
    U = copularnd(f, alpha,N);
    estimate = kendall(U(:,1), U(:,2));
    pass = isnear(estimate, tau, .05);
    if ~pass
        fprintf('For the %s copula, estimated Kendall tau of sample is %.2f while the true value should be %.2f.\n', f,estimate, tau)    
    end
end



    
 