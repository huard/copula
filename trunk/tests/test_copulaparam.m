% TEST COPULAPARAM
fprintf('Test copulaparam ... ')
family = {'Clayton', 'Frank', 'Gumbel', 'Gaussian', 't', 'AMH',  'FGM', 'Arch12', 'Arch14'};
tau = linspace(-.95,.95,10);
for i=1:length(family)
    pass = check_tau(family{i}, tau);
    tau_passed = tau(pass);
    alpha = copulaparam(family{i}, tau_passed);
    ok = isnear(tau_passed, copulastat(family{i}, alpha), 1e-6);
    if any(~ok)
        error('Bug in copulaparam with family ''%s''. Bad [taus; alpha]: \n%s\n%s', ...
            family{i}, num2str(tau_passed(~ok)), num2str(alpha(~ok)) )
    end
end
fprintf('Passed !\n')