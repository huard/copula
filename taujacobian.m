% Function prior_par(type, theta)
% Computes the derivative of function tau = g(theta).
% Type is the copula type.
% theta (parametre) is the vector of the parameters.

function p = prior_par(type,parametre)

if nargin < 2
    error('Requires two input arguments');
end

switch lower(type)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ellipitical copulas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'gaussian'
        rho = parametre;
        
        if (rho < -1 | 1 < rho)
            error('RHO must be a correlation coefficient between -1 and 1');
        else
            p = (2/pi).*(1./sqrt(1-rho.^2));
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Archimedean copulas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case {'clayton' 'frank' 'gumbel' 'fgm' 'amh' 'gb' 'joe' 'arch12' 'arch14'}
        if nargin < 2
            error('Requires two input arguments for an Archimedean copula.');
        end
        alpha = parametre;

        switch lower(type)
            case 'clayton'
                if alpha < -1 and alpha ~= 0
                    error('ALPHA must be nonnegative for the Clayton copula.');
                else
                    p = 2./((alpha + 2).^2);
                end

            case {'frank'} 
                if alpha < 0
                    error('ALPHA must be nonnegative for the Frank copula.');
                else
                    %  t1 = exp(alpha);
                     %  t3 = quadg('dillog',1,t1,1e-1,[]);
                     %   t7 = alpha.^2;
                      %   p = 4 * (-alpha + t1.*alpha - (2 * t3) + 2 * t3 .* t1 + t7 .* t1) ./ (-1 + t1) ./ t7 ./ alpha;
                    load frankdgdata
                    p = interp1(frankdgdata_theta, frankdgdata, abs(alpha), 'spline', 'extrap');
%                    clear frankdgdata frankdgdata_theta
                end



            case 'gumbel'
                if alpha < 1
                    error('ALPHA must be greater than or equal to 1 for the Gumbel copula.');
                else
                    p = 1./(alpha.^2);
                end

            case 'fgm'
                if abs(alpha) > 1
                    error('ALPHA must lie in [-1,1] for the FGM copula.');
                else
                    p = 2/9;
                end

            case 'amh'
                if alpha < -1 or alpha >= 1
                    error('Parameter must be in [-1, 1) for the AMH copula.');
                else
                    t1 = alpha.^ 2;
                    t3 = log(1-alpha);
                    p = -2 / 3 * (t1 + 2 * t3 .* alpha - (2 * alpha) - 2 * t3) ./ t1 ./ alpha;
                end

            case 'arch12'
                if alpha < 1
                    error('ALPHA must be greater than or equal to 1 for the arch12 copula.');
                else
                    p = (2/3)./(alpha.^2);
                end

            case 'arch14'
                if alpha < 1
                    error('ALPHA must be greater than or equal to 1 for the arch14 copula.');
                else
                    p = 4./((1+2*alpha).^2);
                end
        end

    otherwise
        error('Unrecognized copula type: ''%s''',type);
end