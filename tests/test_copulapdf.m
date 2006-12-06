% TEST COPULAPDF
warning off COPULA:BadParameter
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


% Check the calling with 4 arguments
alpha = linspace(-5,5,50);
u = rand(1,10);
v = rand(1,5);
[uu,vv] = meshgrid(u,v);
for i=1:length(families)
    pass = check_alpha(families{i}, alpha);
    alphas = alpha(pass);
    pdf = copulapdf(families{i}, u, v, alphas(4));
    pdf = copulapdf(families{i}, uu, vv, alphas(4));
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