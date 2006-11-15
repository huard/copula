function p = prior_par(type,par)
% Function prior_par(type, theta)
% Computes the derivative of function tau = g(theta) (g corresponds to the function copulastat).
% Type is the copula type. One of {'gaussian' 'clayton' 'frank' 'gumbel' 'fgm' 'amh'
% 'gb' 'joe' 'arch12' 'arch14'}
% par is the vector of the parameters.

if nargin < 2
    error('Requires two input arguments');
end
% Check if the parameters are correct for this copula family
if sum(~bool_par(par,type))>=1;
    [tau_min,tau_max] = int_tau(type);
    error(['PAR must lie in [',num2str(copulaparam(type,tau_min)),',',num2str(copulaparam(type,tau_max)),'] for the ',type,' copula']);
end

switch lower(type)
    case 'gaussian'
        rho = par;
        p = (2/pi).*(1./sqrt(1-rho.^2));
    case {'clayton' 'frank' 'gumbel' 'fgm' 'amh' 'arch12' 'arch14'}
        alpha = par;

        switch lower(type)
            case 'clayton'
                p = 2./((alpha + 2).^2);
            case {'frank'}
                load frankdgdata
                p = interp1(frankdgdata_theta, frankdgdata, abs(alpha), 'spline', 'extrap');
            case 'gumbel'
                p = 1./(alpha.^2);
            case 'fgm'
                p = 2/9.*ones(size(par));
            case 'amh'
                t1 = alpha.^ 2;
                t3 = log(1-alpha);
                p = -2 / 3 * (t1 + 2 * t3 .* alpha - (2 * alpha) - 2 * t3) ./ t1 ./ alpha;
            case 'arch12'
                p = (2/3)./(alpha.^2);
            case 'arch14'
                p = 4./((1+2*alpha).^2);
        end

    otherwise
        error('Unrecognized copula type: ''%s''',type);
end