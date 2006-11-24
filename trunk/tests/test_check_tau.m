% TEST CHECK_TAU
fprintf('Test check_tau ... ')
families = {'amh' 'arch12' 'arch14' 'clayton' 'frank' 'gaussian' 't' 'fgm' 'gumbel'};
for i=1:length(families)
    pass = check_tau(families{i}, [-2,3]);
    if any(pass)
        error('Bug found in check_tau for ''%s''.', families{i})
    end
    pass = check_tau(families{i}, 1/3);
    if any(~pass) & strcmp(families{i}, 'amh')
        error('Bug found in check_tau for ''%s''.', families{i})
    end
end
fprintf('Passed !\n')