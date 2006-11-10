% Données en entrées
load volume.txt
load debit.txt
deb=debit;
vol=volume;
par_deb=[1558.6,460.65];
par_v=[4456.52,915.71];

    % calcul des fractiles
    udeb= normcdf(deb,par_deb(1),par_deb(2));
    uvol = normcdf(vol,par_v(1),par_v(2));
% Octave
%    udeb= normal_cdf(deb,par_deb(1),par_deb(2));
%    uvol = normal_cdf(vol,par_v(1),par_v(2));
%
    u = [udeb';uvol']';

    % calcul des poids associés à chaque famille
    COPULES={'amh' 'gumbel' 'frank' 'clayton' 'arch12' 'arch14' 'gaussian' 'fgm'};
    p = bcs(COPULES, u, [-.97 .97]);
    p = p./sum(p);
    figure;
    bar(real(p));
    set(gca,'XTickLabel',COPULES);
    title(['Poids']);
   
