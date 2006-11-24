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