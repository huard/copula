function c = copula(u1,u2,type,par)
% Function C = copula(u,v,type,par)
%
% Return  C(u,v|par) for the given copula family.
% 
% Input:
% u,v: quantiles
% Type: copula family, one of 'Ind', 'Gumbel', 'Clayton', 'FGM', 'AMH', 'GB', 'Frank', 'Joe'.
% par: copula parameter
% 
% G. Evin, 2005

if (nargin < 3)
    error('Function copula takes a least three arguments.');
end

switch lower(type)
    case 'ind' % Independent
        c = u1.*u2;
        
    case 'gumbel' % Gumbel
        
        if sum(par < 1)>=1
            error('theta >= 1 for Gumbel copula.');
        end
        
        a1 = (-log(u1)).^par;
        a2 = (-log(u2)).^par;
        a3 = 1./par;
        c = exp(-(a1+a2).^a3);
        
    case 'clayton'
        if sum(par < 0)>=1
            error('theta >= 0 for Clayton copula.');
        end
        c = (u1.^(-par)+u2.^(-par)-1).^(-1./par);
        
    case 'fgm' % Farlie, Gumbel, Morgenstern
        c = u1.*u2.*(1+par.*(1-u1).*(u2));
        
    case 'amh' % Ali, Mikhail, Haq
        if sum(par < -1 || par >= 1)>=1
            error('theta in [-1,1[ for Ali, Mikhail and Haq copula.');
        end
        c = (u1.*u2)./(1-par.*(1-u1).*(u2));
        
    case 'gb' % Gumbel & Barnett
        if sum(par < 0 || par > 1)>=1
            error('theta in [0,1] for Gumbel & Barnett copula.');
        end
        c = u1.*u2.*exp(-par.*log(u1).*log(u1));
        
    case 'frank' % Frank, 1979
        if sum(par == 0)>=1
            error('theta != 0 for Frank copula.');
        end
        c = (-1./par).*log(1+(exp(-par.*u1)-1).*(exp(-par.*u2)-1)./(exp(-par)-1));
        
    case 'joe' % Joe, 1997
        if sum(par < 1)>=1
            error('theta >= 1 for Joe copula.');
        end
        c = 1 - (((1-u1).^par)+((1-u2).^par)-((1-u1).^par).*((1-u2).^par)).^(1./par);
        
    otherwise
        error('Unknown copula family:',type);
end
