% Test copularnd

family = {'ind', 'gaussian', 't', 'clayton', 'frank', 'gumbel',... 
		        'fgm', 'amh', 'arch12', 'arch14', 'gb', 'joe'}

% Test correlation of generated sample.
% It should match the tau corresponding to the parameter given to
% copularnd.
% N = 500
% for i=1:length(family)
%     f = family{i};
%     bound= tauboundaries(f);
%     tau = bound(1) + rand*diff(bounds);
%     alpha = copulaparam(f, tau);
%     U = copularnd(f, alpha,N);
%     pass(i) = isnear(kendall(U), tau, 1e-2);
% end
% if any(~pass)
%     warning('Copularnd seem to have a problem with some families.\n')
%     fprintf('%s ', family(~pass))
% end


% Overlay a random sample with copula density.
family = {'gaussian', 'gumbel' 'clayton' 'frank' 'amh' 'fgm' 'arch12' 'arch14'};     
N = 3000;
M = 100;
for i=1:length(family)
    f = family{i};
    bounds = tauboundaries(f);
    tau = bounds(1) + .7*diff(bounds);
    alpha = copulaparam(f, tau);
    U = copularnd(f, alpha,N);
    x = linspace(0,1,M);
    [u,v] = meshgrid(x,x);
    u1 = reshape(u, M^2, 1);
    v1 = reshape(v, M^2, 1);
    pdf = copulapdf(f, [u1,v1], alpha);
    pdf = reshape(pdf, M, M);
    subplot(3,3,i)
    axis equal
    contour(u,v,pdf,50);
    hold on
    scatter(U(:,1), U(:,2), 2)
    title(f)
    hold off
end
saveas(gcf, 'density', 'png')
close
    

    
 