function y = generateurarch(x,type,parametre)
%function y = generateurarch(x,type,parametre)
%
% Function generator of the archimedean copulas
%
% INPUTS:   x is a vector that lies in [0,1]
%           type is one of the list {'ind' 'gumbel' 'clayton' 'frank' 'gb' 'amh' 'joe'}
%           parametre is the parameter of the copula, none for the independant copula
%
% OUTPUTS:	y, a vectors that contains the results of the generator on the
% points x
%
% Guillaume EVIN
%
%  13 May, 2004.

if (nargin < 2 || nargin > 3)
    error('la fonction nécessite deux ou trois arguments');
end

y = ones(size(x)).*realmax; %le générateur peut tendre vers l'infini en 0
nonnul = (x ~= 0);
    
switch lower(type)
    
    case 'ind' %independance
        y(nonnul) = -log(x(nonnul));
        
    case 'amh' %Copule de Ali, Mikhail and Haq 1978
        if (parametre < -1 || parametre >= 1)
            error('on doit avoir theta dans [-1,1[ pour le copule de Ali, Mikhail et Haq');
        end
        y(nonnul) = log((1-parametre*(1-x(nonnul)))./x(nonnul));
        
    case 'gb' %Copule de Gumbel 1960, Barnett 1980
        if (parametre < 0 || parametre > 1)
            error('on doit avoir theta dans [0,1] pour le copule de Gumbel et Barnett');
        end
        y(nonnul) = log(1-parametre*log(x(nonnul))); 
        
        
    case 'gumbel' %Copule de Gumbel 1960, Hougaard 1986
        if parametre < 1
            error('on doit avoir theta >= 1 pour le copule de Gumbel');
        end
        y(nonnul) = (-log(x(nonnul))).^parametre; 
        
    case 'clayton' %Copule de Clayton (Clayton 1978)
        if parametre < 0
            error('on doit avoir theta >= 0 pour le copule de Clayton');
        end
        y(nonnul) = x(nonnul).^(-parametre)-1;
        
    case 'frank' %Copule de Frank (Frank 1979)
        if parametre == 0
            error('on doit avoir theta different de 0 pour le copule de Franck');
        end
        y(nonnul) = -log((exp(-parametre.*x(nonnul))-1)./(exp(-parametre)-1)); 
        
    case 'joe' %Copule de Joe (Joe 1997)
        if parametre < 1
            error('on doit avoir theta >= 1 pour le copule de Joe');
        end
        
        y(nonnul) = -log(1-(1-x(nonnul)).^parametre);  
        
    
otherwise
    error('Type de copule archimédien inconnu: ',type);
end