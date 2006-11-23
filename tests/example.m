% Example of application in hydrology
% Finding the copula family that describes best the correlation between
% maximum annual river discharges and volumes. 

% Load data
load volume.txt
load discharge.txt

% Compute quantiles assuming normal margins. 
par_dis=[1558.6,460.65];
par_vol=[4456.52,915.71];
u= normcdf(discharge,par_dis(1),par_dis(2));
v = normcdf(volume,par_vol(1),par_vol(2));

% Data
U = [u';v']';

% Prior on tau (a little awkward to maintain compatibility with MatLab 6).
prior_tau = inline('1/sqrt(2*pi)/.3 * exp(-((x-.5).^2)/2/.3^2)', 'x');

% Copula families to test.
COPULAS={'AMH' 'Gumbel' 'Frank' 'Clayton' 'Arch12' 'Arch14' 'Gaussian' 'FGM' 'Ind'};

% Compute weight of each copula family.
p = bcs(COPULAS, U, [-.97 .97], prior_tau);

% This would be the place to include the effect of the prior on the copula
% family, if different from the one defined in the paper by Huard et al.

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

