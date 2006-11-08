function pass = check_domain_tau(family, theta)

% Function pass = check_domain_tau(family, theta)
%
% Return boolean array 
% True if theta is in the family domain
% False if not.
%

% D. Huard, Nov. 2006

tau = copulastat(family, theta)
bounds = tauboundaries(family)
n = size(bounds)
pass = zeros(size(theta))
for i=1:n(1)
    pass =  pass | ((tau > bounds(i,1)) & (tau < bounds(i,2)))
end


