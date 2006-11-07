function y = generateurarchderivee2(x,type,parametre)
%function y = generateurarchderivee(x,type,parametre)
%
% Function second derivate of the generator of the archimedean copulas, positive on [0,1]
%
% INPUTS:   x is a vector that lies in [0,1]
%           type is one of the list {'ind' 'gumbel' 'clayton' 'frank' 'gb' 'amh' 'joe'}
%           parametre is the parameter of the copula, none for the independant copula
%
% OUTPUTS:	y, a vectors that contains the second derivate of the generator on the
% points x
%
% Guillaume EVIN
%
%  13 May, 2004.

if (nargin < 2 || nargin > 3)
    error('la fonction nécessite deux ou trois arguments');
end

if sum(x < 0) > 0 || sum(x > 1) > 0
 error('un generateur archimedien n est defini que sur [0,1]');
end

y = ones(size(x)).*realmax;
nonnul = (x ~= 0);

switch lower(type)
    
    case 'ind' %independance
        y(nonnul) = 1./((x(nonnul).^2));
        
    case 'amh' %Copule de Ali, Mikhail and Haq 1978
        if (parametre < -1 || parametre >= 1)
            error('on doit avoir theta dans [-1,1[ pour le copule de Ali, Mikhail et Haq');
        end
        y(nonnul) = (-parametre^2)./((1-parametre.*(1-x(nonnul))).^2)+1./(x(nonnul).^2);
        
    case 'gb' %Copule de Gumbel 1960, Barnett 1980
        if (parametre < 0 || parametre > 1)
            error('on doit avoir theta dans [0,1] pour le copule de Gumbel et Barnett');
        end
        f = 1 - parametre.*log(x(nonnul));
        y(nonnul) = (parametre.*(f-parametre))./((x(nonnul).^2).*(f.^2)); 
        
    case 'gumbel' %Copule de Gumbel 1960, Hougaard 1986
        if parametre < 1
            error('on doit avoir theta >= 1 pour le copule de Gumbel');
        end
        y(nonnul) = parametre.*parametre.*((-log(x(nonnul))).^(parametre-2))./(x(nonnul).^2); 
        
    case 'clayton' %Copule de Clayton (Clayton 1978)
        if parametre < 0
            error('on doit avoir theta >= 0 pour le copule de Clayton');
        end
        y(nonnul) = (parametre*(parametre+1)).*(x(nonnul).^(-parametre-2));
        
    case 'frank' %Copule de Frank (Frank 1979)
        if parametre == 0
            error('on doit avoir theta different de 0 pour le copule de Franck');
        end
        f = exp(-parametre.*x(nonnul));
        y(nonnul) = ((parametre^2).*f)./((f-1).^2);
        
    case 'joe' %Copule de Joe (Joe 1997)
        if parametre < 1
            error('on doit avoir theta >= 1 pour le copule de Joe');
        end
        f = (1 - x(nonnul));
        y(nonnul) = (parametre.*(f.^(parametre-2)).*(f.^parametre-1+parametre))./((f.^parametre-1).^2);  
        
        
    otherwise
        error('Type de copule archimédien inconnu: ',type);
end