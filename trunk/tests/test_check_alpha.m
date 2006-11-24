% TEST CHECK_ALPHA
fprintf('Test check_alpha ... ')
family = {'amh' 'arch12' 'arch14' 'clayton' 'frank' 'gaussian' 'fgm' 'gumbel'};
for i=1:length(family)
    bounds = tauboundaries(family{i});
    tau = linspace(bounds(1)+.01, bounds(2)-.01, 10);
    alpha = copulaparam(family{i}, tau);
    pass = check_alpha(family{i}, alpha);
    if any(~pass)
        error('Bug in check_alpha for ''%s''.', family{i})
    end
end
if ~check_alpha('t', -.9, 2) | check_alpha('t', -.9, 2.4)
    error('Bug in check_alpha for t.')
end
try
 check_alpha('t', -.9);
 catch
 end
fprintf('Passed !\n')
