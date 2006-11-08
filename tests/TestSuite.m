% TEST SUITE

% TEST CHECK_TAU
families = {'amh' 'arch12' 'arch14' 'clayton' 'frank' 'gaussian' 't' 'fgm' 'gumbel'};
for i=1:length(families)
    pass = check_tau(families{i}, [-2,3],0);
    if any(pass)
        error('Bug found in check_tau for ''%s''.', families{i})
    end
    pass = check_tau(families{i}, 1/3, 0);
    if any(~pass) & families{i}~='amh'
        error('Bug found in check_tau for ''%s''.', families{i})
    end
end


% TEST CHECK_ALPHA
% TEST COPULASTAT
families = {'Clayton', 'Gumbel', 'Gaussian', 't', 'AMH', 'GB', 'Joe', 'FGM', 'Arch12', 'Arch14'};
for i=1:length(families)
    try 
        tau = copulastat(families{i}, -2);
        error('This should raise an error: copulastat(%s, -2)', families{i})
    catch
        %ok
    end
end
   