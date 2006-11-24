% TEST COPULASTAT
fprintf('Test copulastat ... ')
families = {'Clayton', 'Gumbel', 'Gaussian', 'AMH', 'FGM', 'Arch12', 'Arch14', 'Frank'};
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

%% t
if any(~isnear([-1,1.047197551/pi, 2.239539030/pi], copulastat('gaussian', [-1, .5, .9]), 1e-6))
    error('Bug in copulastat for t.')
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
fprintf('Passed !\n')
