% TEST SUITE FOR COPULA RELATED FUNCTIONS
%
% Modified by:
%% D. Huard,  Nov. 8, 2006. Added TEST CHECK_TAU, TEST CHECK_ALPHA, TEST
%% COPULASTA

% TEST CHECK_TAU
fprintf('Test check_tau ... ')
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
fprintf('Passed !\n')



% TEST CHECK_ALPHA
fprintf('Test check_alpha ... ')
alpha = [1,2,3,1000];
pass = check_alpha('gumbel', alpha);
if any(~pass)
    error('Bug in check_alpha for Gumbel')
end
fprintf('Passed !\n')

warning off COPULA:BadParameter
% TEST COPULASTAT
fprintf('Test copulastat ... ')
families = {'Clayton', 'Gumbel', 'Gaussian', 't', 'AMH', 'GB', 'Joe', 'FGM', 'Arch12', 'Arch14', 'Frank'};
alpha = linspace(-5,5,50);
for i=1:length(families)
    pass = check_alpha(families{i}, alpha);
    tau = copulastat(families{i}, alpha(pass));
end
for i=1:length(families)
    try
        tau = copulastat(families{i}, -2)
        fprintf('This should raise an error: copulastat(%s, -2)', families{i})
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
% compare copulastat('gumbel', alpha);
fprintf('Passed !\n')


% TEST COPULAPARAM

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
fprintf('Passed !\n')

% TEST COPULACDF
% Values are compared with results obtained with Maple.

%% Clayton
pass = isnear([0.206897, 0.247226, 0.268499], copulacdf('clayton', [.3,.4], [1,2,3]), 1e-6)
if any(~pass):
    error('Error in Clayton copulacdf.')
end

%% Gumbel
pass = isnear([.12, .22025, .256704, .297721], copulacdf('gumbel', [.3,.4], [1,2,3,10]), 1e-6)
if any(~pass):
    error('Error in Gumbel copulacdf.')
end

%% Frank
pass = isnear([0.4539769530e-2, .1452283617, .1693783396, 0.1911490534], copulacdf('frank', [.3,.4], [-10, 1, 2,3]), 1e-9)
if any(~pass):
    error('Error in Frank copulacdf.')
end



% TEST COPULA_LIKE

% TEST TAUJACOBIAN
fprintf('Test taujacobian ... ')
families = {'Clayton', 'Gumbel', 'Gaussian', 'AMH', 'FGM', 'Arch12', 'Arch14', 'Frank'};
for i=1:length(families)
    pass = check_alpha(families{i}, alpha);
    j = taujacobian(families{i}, alpha(pass));
end
fprintf('Passed !\n')

% TEST BCS


