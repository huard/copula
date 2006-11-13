% TEST SUITE FOR COPULA RELATED FUNCTIONS
% Some tests check that the functions return matrices of the correct shape.
% Other functions check that the actual values are vaild, by comparing
% results with those obtained using Maple. See file analytic.mws to see how
% these computations were performed.

% Modified by:
%% D. Huard,  Nov. 8, 2006. Added TEST CHECK_TAU, TEST CHECK_ALPHA
%% D. Huard, Nov, 9, 2006. Added TEST COPULAPDF, TEST COPULACDF, COPULA
%%  TEST COPULAPARAM, TEST COPULASTAT
%%  Found bug in copulacdf('clayton')
%% D. Huard, Nov, 12, 2006. Added TAUJACOBIAN
%%  Found bug in copulacdf('fgm')

warning('off');%, COPULA:BadParameter)

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



% TEST CHECK_ALPHA
fprintf('Test check_alpha ... ')
alpha = [1,2,3,1000];
pass = check_alpha('gumbel', alpha);
if any(~pass)
    error('Bug in check_alpha for Gumbel')
end
fprintf('Passed !\n')


% TEST COPULASTAT
fprintf('Test copulastat ... ')
families = {'Clayton', 'Gumbel', 'Gaussian', 't', 'AMH', 'FGM', 'Arch12', 'Arch14', 'Frank', 'GB', 'Joe'};
alpha = linspace(-5,5,50);
for i=1:length(families)
    pass = check_alpha(families{i}, alpha);
    tau = copulastat(families{i}, alpha(pass));
end
for i=1:length(families)
    try
        tau = copulastat(families{i}, -2);
        if ~strcmp(families{i}, 'Frank')
            fprintf('This should raise an error: copulastat(%s, -2)', families{i})
        end
    catch
    end
end
try
    tau = copulastat('frank', 0)
    fprintf('This should raise an error: copulastat(''frank'', 0)')
catch 
end
try 
    tau = copulastat('frankie', .5)
    fprintf('This should raise an error: copulastat(''frankie'', 0)')
catch
end

%% Clayton
if any(~isnear([0, .5, 3/4], copulastat('clayton', [0,2,6]), 1e-6))
    error('Bug in copulastat for Clayton.')
end

%% Gumbel
if any(~isnear([0,.5, .8], copulastat('gumbel', [1,2,5]), 1e-6))
    error('Bug in copulastat for Gumbel.')
end

%% Frank
if any(~isnear([-.6657773860, .1100185380, .2138945700, .4567009584], copulastat('frank', [-10, 1,2,5]), 1e-6))
    error('Bug in copulastat for Frank.')
end

%% AMH
if any(~isnear([5/3 - 8/3*log(2),-0.099457315, 0.3269125715], copulastat('amh', [-1,-.5,.99]), 1e-6))
    error('Bug in copulastat for AMH.')
end

%% Gaussian
if any(~isnear([-1,1.047197551/pi, 2.239539030/pi], copulastat('gaussian', [-1, .5, .9]), 1e-6))
    error('Bug in copulastat for Gaussian.')
end

%% Arch12
if any(~isnear([1/3, 2/3, 8/9], copulastat('arch12', [1,2,6]), 1e-6))
    error('Bug in copulastat for Arch12.')
end

%% Arch14
if any(~isnear([1/3, 1/2, 4/5], copulastat('arch14', [1, 1.5, 4.5]), 1e-6))
    error('Bug in copulastat for Arch14.')
end

%% FGM
if any(~isnear([-2/9, 0, 2/27, 2/9], copulastat('fgm', [-1, 0, 1/3, 1]), 1e-6))
    error('Bug in copulastat for FGM.')
end
fprintf('Passed !')
fprintf('  Test missing for Joe and GB. \n')

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


% TEST COPULAPDF
fprintf('Test copulapdf ... ')
families = {'ind', 'gaussian', 'gumbel' 'clayton' 'frank' 'amh' 'joe' 'fgm' 'arch12' 'arch14'};
pdf = copulapdf('gumbel', [[.3,.4];[.5,.6]], [4,5,6]);
if ~all(size(pdf) == [2,3])
    error('Bad shape for array returned by copulapdf.')
end

pdf = copulapdf('clayton', [[.3,.4]; [.5,.6]], [4;5]);
if ~all(size(pdf) == [2, 1])
    error('Bad shape for array returned by copulapdf.')
end

try 
    pdf = copulapdf('gumbel', [[.3,.4],[.5,.6]], [4;5;6]);
    fprintf('This should raise an error due to bad arguments.')
catch
end
   
alpha = [[4,5,6];[-1,-2,-3]];
pdf = copulapdf('frank',[[.8,.9]], alpha);
if size(pdf) ~= size(alpha)
    error('Bad shape for array returned by copulapdf.')
end

alpha = linspace(-5,5,50);
U = [[.3,.4];[.5,.6]];
for i=1:length(families)
    pass = check_alpha(families{i}, alpha);
    pdf = copulapdf(families{i}, U, alpha(pass));
end

% Check analytical values with Maple
%% Clayton
pass = isnear([1.003064594, 1.603413484, 2.237830668], copulapdf('clayton', [.3,.4], [.1,2,5]), 1e-8);
if any(~pass)
    error('Error in Clayton copulapdf.')
end

%% Gumbel
pass = isnear([1., 1.469156047, 1.881228154, 1.598617142], copulapdf('gumbel', [.3,.4], [1,2,3,10]), 1e-8);
if any(~pass)
    error('Error in Gumbel copulapdf.')
end

%% Frank
pass = isnear([.4546784123, 1.050358851, 1.450640692, 0.4539580774e-2], copulapdf('frank', [.3,.4], [-10, 1, 5,100]), 1e-8);
if any(~pass)
    error('Error in Frank copulapdf.')
end

%% AMH
pass = isnear([1.007941249, 1.024223647, 1.156074898], copulapdf('amh', [.3,.4], [.1, .3, .9]), 1e-9);
if any(~pass)
    error('Error in AMH copulapdf.')
end

%% Gaussian
pass = isnear([0.5933685602, 1.192296359, 1.240318192], copulapdf('gaussian', [.3,.4], [-.9,.5,.99]), 1e-9);
if any(~pass)
    error('Error in Gaussian copulapdf.')
end

%% Arch12
pass = isnear([1.230062734, 1.992107671, 0.2548015110], copulapdf('arch12', [.3,.4], [1,5,12]), 1e-8);
if any(~pass)
    error('Error in Arch12 copulapdf.')
end
                          
%% Arch14
pass = isnear([1.230062733, 2.321761087, 1.061073261], copulapdf('arch14', [.3,.4], [1,5,12]), 1e-7);
if any(~pass)
    error('Error in Arch14 copulapdf.')
end                   
                          
%% FGM
pass = isnear([.928, 1.008, 1.024, 1.072], copulapdf('fgm', [.3,.4], [-.9, .1, .3, .9]), 1e-9);
if any(~pass)
    error('Error in FGM copulapdf.')
end

fprintf('Passed !\n')

% TEST COPULACDF
% Values are compared with results obtained with Maple.
fprintf('Test copulacdf ... ')

%% Clayton
pass = isnear([0.1325957657, 0.2472256930, 0.2876051623], copulacdf('clayton', [.3,.4], [.1,2,5]), 1e-8);
if any(~pass)
    error('Error in Clayton copulacdf.')
end

%% Gumbel
pass = isnear([.12, .22025, .256704, .297721], copulacdf('gumbel', [.3,.4], [1,2,3,10]), 1e-6);
if any(~pass)
    error('Error in Gumbel copulacdf.')
end

%% Frank
pass = isnear([0.4539769530e-2, .1452283617, .2255806654, .2995464418], copulacdf('frank', [.3,.4], [-10, 1, 5, 40]), 1e-6);
if any(~pass)
    error('Error in Frank copulacdf.')
end

%% AMH
pass = isnear([.1252609603, 0.1372997712, .1929260450 ], copulacdf('amh', [.3,.4], [.1, .3, .9]), 1e-9);
if any(~pass)
    error('Error in AMH copulacdf.')
end

%% Arch12
pass = isnear([.2068965517, .2956430471, .2999130447], copulacdf('arch12', [.3,.4], [1,5,12]), 1e-8);
if any(~pass)
    error('Error in Arch12 copulacdf.')
end

%% Arch14
pass = isnear([.2068965517, .2873095536, .2990807552], copulacdf('arch14', [.3,.4], [1,5,12]), 1e-7);
if any(~pass)
    error('Error in Arch14 copulacdf.')
end    

%% FGM
pass = isnear([.07464, .12504,.13512, .16536], copulacdf('fgm', [.3,.4], [-.9, .1, .3, .9]), 1e-8);
if any(~pass)
    error('Error in FGM copulacdf.')
end

fprintf('Passed !\n')


% TEST COPULA_LIKE

% TEST TAUJACOBIAN
fprintf('Test taujacobian ... ')
families = {'Clayton', 'Gumbel', 'Gaussian', 'AMH', 'FGM', 'Arch12', 'Arch14', 'Frank'};
for i=1:length(families)
    pass = check_alpha(families{i}, alpha);
    j = taujacobian(families{i}, alpha(pass));
end


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

% TEST LAMBDAARCH
family = {'Clayton', 'Frank', 'Gumbel', 'AMH', 'Joe', 'GB'};
alpha = linspace(-5,5,10);
for i=1:length(family)
    pass = check_alpha(family{i}, alpha);
    alpha_passed = alpha(pass);
    taus1 = copulastat(family{i}, alpha_passed);
    taus2 = zeros(size(taus1));
    for j=1:length(taus1)
        try
            taus2(j) = 1 + 4 .* quadg('lambdaarch',0,1,[],[],family{i}, alpha_passed(j));
        catch
            family{i}
        end
    end
    if any(~isnear(taus1, taus2, 1e-4))
        error('Bug in lambdaarch for copula ''%s''.\ncopulastat: %s\nlambdaarch: %s', family{i}, num2str(taus1), num2str(taus2))
    end
end 

fprintf('Everything looks fine.\n')

% TEST BCS
run TEST_acf
%close


