% TEST TAUJACOBIAN
% First check that p(alpha) is normalized properly. 
warning off MATLAB:fzero:UndeterminedSyntax

fprintf('Test taujacobian ... ')
families = {'Clayton', 'Gumbel', 'Gaussian', 'AMH', 'FGM', 'Arch12', 'Arch14', 'Frank'};
for i=1:length(families)
    family = families{i};
    taus = tauboundaries(family);
    alphas = copulaparam(family, [taus(1)+.001, taus(2)-.006]);
    integral = quadg('tau_jacobian2',alphas(1), alphas(2), 1e-5, [0,256], family)/diff(taus);
    if integral < .98
        error('The prior on alpha does not properly normalized for copula ''%s''.', family)
    end
end

% Check values are similar to analytically computed values.
if any(~isnear([1/2, 2/9, 1/8], taujacobian('clayton', [0, 1, 2]), 1e-8))
    error('Bug in taujacobian for Clayton.')
end

if any(~isnear([1,1/4, 1/9], taujacobian('gumbel', [1,2,3]), 1e-8))
    error('Bug in taujacobian for Gumbel.')
end

if any(~isnear([0.2686268354e-1, .1078697551, 0.6274654072e-1, 0.3868405077e-3], taujacobian('frank', [-10,1,5,100]), 1e-9))
    error('Bug in taujacobian for Frank.')
end

if any(~isnear([2.010075630/pi, 2.309401076/pi, 4.588314676/pi], taujacobian('gaussian', [-.1, .5, .9]), 1e-6))
    error('Bug in taujacobian for Gaussian.')
end

if any(~isnear([2-8/3*log(2), .3032150400, .4842094029], taujacobian('amh', [-1, .5, .9]), 1e-6))
    error('Bug in taujacobian for AMH.')
end

if any(~isnear([2/3, 1/6, 2/27], taujacobian('arch12', [1,2,3]), 1e-6))
    error('Bug in taujacobian for Arch12.')
end

if any(~isnear([4/9, 4/25, 4/49], taujacobian('arch14', [1,2,3]), 1e-6))
    error('Bug in taujacobian for Arch14.')
end

if any(~isnear([2/9,2/9,2/9,2/9], taujacobian('fgm', [-1,0,.5,1]), 1e-6))
    error('Bug in taujacobian for FGM.')
end
fprintf('Passed !\n')