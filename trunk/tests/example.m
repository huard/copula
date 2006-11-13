% Example
% Finding the copula family that describes best the correlation between
% maximum annual discharges and volumes. 

% Load data
load volume.txt
load discharge.txt

% Compute fractile assuming normal margins. 
par_dis=[1558.6,460.65];
par_vol=[4456.52,915.71];

u= normcdf(discharge,par_dis(1),par_dis(2));
v = normcdf(volume,par_vol(1),par_vol(2));

% Octave
%    udeb= normal_cdf(deb,par_deb(1),par_deb(2));
%    uvol = normal_cdf(vol,par_v(1),par_v(2));

U = [u';v']';

% Compute weight of each copula family.
COPULAS={'AMH' 'Gumbel' 'Frank' 'Clayton' 'Arch12' 'Arch14' 'Gaussian' 'FGM'};
p = bcs(COPULAS, U, [-.97 .97]);

% This is the place where you add the effect of the prior on the copula
% family, if different from uniform.

% Normalize weights
p = p./sum(p);

% Print results
fprintf('\n%-10s\t%-10s', 'Family', 'Weight')
for i=1:length(COPULAS)
    fprintf('\n%-10s\t%-10f', COPULAS{i}, p(i))
    if p(i) == max(p)
        fprintf(' <==')
    end
end

% Plot results
figure;
bar(real(p));
set(gca,'XTickLabel',COPULAS);
title(['Weight']);

