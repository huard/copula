function u = lambdaarch(t, family, alpha)
%
%    FUNCTION U = LAMBDAARCH(T, FAMILY, ALPHA)
%
%    Function used in the paper of Genest [1993], in order to compare
%    differents copulas, lambda(t) = g(t)/g'(t) where g is the generator
%
%    INPUTS   
%       t:      Vector that lies in [0,1]
%       FAMILY: One of {'gumbel' 'clayton' 'frank' 'gb' 'amh' 'joe'} 	
%       ALPHA:  Parameter of the copula
%
%    OUTPUTS	
%       U:      Vector of lambda at points t.
%
% Guillaume EVIN
%
%  13 May, 2004.

if strcmp(lower(family), 'ind')
    if nargin ~= 2
        error('Function takes two arguments.')
   end
else
    if nargin ~=3
        error('Function takes three arguments.')
    end
end

pass = check_alpha(family, alpha);
if any(~pass)
    error('Bad parameters: %s', num2str(alpha(~pass)))
end

den = archgender(t,family,alpha);
nz = (den~=0);
z = (den == 0);

u(nz) = archgen(t(nz),family,alpha)./archgender(t(nz),family,alpha);
u(z) = sign(archgen(t(z),family,alpha)).*realmax;

function y = archemedeangenerator2(family, x, alpha)
%
%    function y = archemedeangenerator(family, x[, alpha])
%
%    Function generator of the archimedean copulas
%
%    INPUT   
%	X:      Vector that lies in [0,1]
%       FAMILY: One of  {'ind' 'gumbel' 'clayton' 'frank' 'gb' 'amh' 'joe'}
%       ALPHA:  Parameter of the copula (None for independent).
%
%    OUTPUT
%	Y:      Vector of the generator at points x.
%
% Guillaume EVIN
%
%  13 May, 2004.



% Generator may tend to infinity at 0
y = ones(size(x)).*realmax; 
nonzero = (x ~= 0);

switch lower(family)
    
    case 'ind' % independent
        y(nonzero) = -log(x(nonzero));
        
    case 'amh' % Ali, Mikhail and Haq 1978
        if (alpha < -1 || alpha >= 1)
            error('on doit avoir theta dans [-1,1[ pour le copule de Ali, Mikhail et Haq');
        end
        y(nonzero) = log((1-alpha*(1-x(nonzero)))./x(nonzero));
        
    case 'gb' % Gumbel 1960, Barnett 1980
        y(nonzero) = log(1-alpha*log(x(nonzero))); 
        
        
    case 'gumbel' %Copule de Gumbel 1960, Hougaard 1986
        y(nonzero) = (-log(x(nonzero))).^alpha; 
        
    case 'clayton' %Copule de Clayton (Clayton 1978)
        y(nonzero) = x(nonzero).^(-alpha)-1;
        
    case 'frank' %Copule de Frank (Frank 1979)
        y(nonzero) = -log((exp(-alpha.*x(nonzero))-1)./(exp(-alpha)-1)); 
        
    case 'joe' %Copule de Joe (Joe 1997)       
        y(nonzero) = -log(1-(1-x(nonzero)).^alpha);  
        
    otherwise
        error('Unrecognized copula family: ''%s''',family);
end

function y = archemedeangeneratorderivative2(family, x, alpha);
%
%    function y = archemedeangeneratorderivative(family, x, alpha)
%
%    Function second derivate of the generator of the archimedean copulas, positive on [0,1]
%
%   INPUT   
%	    FAMILY: One of  {'ind' 'gumbel' 'clayton' 'frank' 'gb' 'amh' 'joe'}
%	    X:      Vector that lies in [0,1]
%       ALPHA:  Parameter of the copula (None for independent).
%
%    OUTPUT
%	    Y:      Vector of the second derivative of generator at points x.
%

% Guillaume EVIN
%  13 May, 2004.

y = ones(size(x)).*realmax;
nonzero = (x ~= 0);

switch lower(family)
    
    case 'ind' % independent
        y(nonzero) = 1./((x(nonzero).^2));
        
    case 'amh' % Ali, Mikhail and Haq 1978
           y(nonzero) = (-alpha.^2)./((1-alpha.*(1-x(nonzero))).^2) + 1./(x(nonzero).^2);
        
    case 'gb' %Gumbel 1960, Barnett 1980
        f = 1 - alpha.*log(x(nonzero));
        y(nonzero) = (alpha.*(f-alpha))./((x(nonzero).^2).*(f.^2)); 
        
    case 'gumbel' % Gumbel 1960, Hougaard 1986
        y(nonzero) = alpha.*alpha.*((-log(x(nonzero))).^(alpha-2))./(x(nonzero).^2); 
        
    case 'clayton' % Clayton (Clayton 1978)
        y(nonzero) = (alpha*(alpha+1)).*(x(nonzero).^(-alpha-2));
        
    case 'frank' % Frank (Frank 1979)
        f = exp(-alpha.*x(nonzero));
        y(nonzero) = ((alpha^2).*f)./((f-1).^2);
        
    case 'joe' % Joe (Joe 1997)
        f = (1 - x(nonzero));
        y(nonzero) = (alpha.*(f.^(alpha-2)).*(f.^alpha-1+alpha))./((f.^alpha-1).^2);  
        
    otherwise
        error('Unrecognized copula family: ''%s''',family);
end

function y = archemedeangeneratorderivative(x, family, alpha)
%
%    function y = archemedeangeneratorderivative(family, x, alpha)
%
%    Function second derivate of the generator of the archimedean copulas, positive on [0,1]
%
%   INPUT   
%	    FAMILY: One of  {'ind' 'gumbel' 'clayton' 'frank' 'gb' 'amh' 'joe'}
%	    X:      Vector that lies in [0,1]
%       ALPHA:  Parameter of the copula (None for independent).
%
%    OUTPUT
%	    Y:      Vector of the derivative of generator at points x.
%

% Guillaume EVIN
%  13 May, 2004.

if (nargin < 2 || nargin > 3)
    error('la fonction nécessite deux ou trois arguments');
end

y = ones(size(x)).*-realmax;
nonzero = (x ~= 0);

switch lower(family)
    
    case 'ind' %independance
        y(nonzero) = -1./x(nonzero);
        
    case 'amh' %Copule de Ali, Mikhail and Haq 1978
        if (alpha < -1 || alpha >= 1)
            error('on doit avoir theta dans [-1,1[ pour le copule de Ali, Mikhail et Haq');
        end
        y(nonzero) = alpha./(1-alpha.*(1-x(nonzero)))-1./x(nonzero);
                
    case 'gb' %Copule de Gumbel 1960, Barnett 1980
        if (alpha < 0 || alpha > 1)
            error('on doit avoir theta dans [0,1] pour le copule de Gumbel et Barnett');
        end
        y(nonzero) = (alpha./x(nonzero))./(alpha.*log(x(nonzero))-1); 
        
    case 'gumbel' %Copule de Gumbel 1960, Hougaard 1986
        if alpha < 1
            error('on doit avoir theta >= 1 pour le copule de Gumbel');
        end
        
        y(nonzero) = (-alpha./x(nonzero)).*((-log(x(nonzero))).^(alpha-1)); 
                
    case 'clayton' %Copule de Clayton (Clayton 1978)
        if alpha < 0
            error('on doit avoir theta >= 0 pour le copule de Clayton');
        end
        y(nonzero) = (-alpha).*(x(nonzero).^(-alpha-1));
            
    case 'frank' %Copule de Frank (Frank 1979)
        if alpha == 0
            error('on doit avoir theta different de 0 pour le copule de Franck');
        end
        y(nonzero) = (alpha.*exp(-alpha.*x(nonzero)))./(-1+exp(-alpha.*x(nonzero)));
                
    case 'joe' %Copule de Joe (Joe 1997)
        if alpha < 1
            error('on doit avoir theta >= 1 pour le copule de Joe');
        end
        y(nonzero) = (alpha.*((1-x(nonzero)).^(alpha-1)))./((1-x(nonzero)).^alpha-1);
            
    otherwise
        error('family de copule archimédien inconnu: ',family);
end

function y = archemedeangenerator(x,family,alpha)
%function y = generateurarch(x,family,alpha)
%
% Function generator of the archimedean copulas
%
% INPUTS:   x is a vector that lies in [0,1]
%           family is one of the list {'ind' 'gumbel' 'clayton' 'frank' 'gb' 'amh' 'joe'}
%           alpha is the parameter of the copula, none for the independant copula
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
nonzero = (x ~= 0);
    
switch lower(family)
    
    case 'ind' %independance
        y(nonzero) = -log(x(nonzero));
        
    case 'amh' %Copule de Ali, Mikhail and Haq 1978
        if (alpha < -1 || alpha >= 1)
            error('on doit avoir theta dans [-1,1[ pour le copule de Ali, Mikhail et Haq');
        end
        y(nonzero) = log((1-alpha*(1-x(nonzero)))./x(nonzero));
        
    case 'gb' %Copule de Gumbel 1960, Barnett 1980
        if (alpha < 0 || alpha > 1)
            error('on doit avoir theta dans [0,1] pour le copule de Gumbel et Barnett');
        end
        y(nonzero) = log(1-alpha*log(x(nonzero))); 
        
        
    case 'gumbel' %Copule de Gumbel 1960, Hougaard 1986
        if alpha < 1
            error('on doit avoir theta >= 1 pour le copule de Gumbel');
        end
        y(nonzero) = (-log(x(nonzero))).^alpha; 
        
    case 'clayton' %Copule de Clayton (Clayton 1978)
        if alpha < 0
            error('on doit avoir theta >= 0 pour le copule de Clayton');
        end
        y(nonzero) = x(nonzero).^(-alpha)-1;
        
    case 'frank' %Copule de Frank (Frank 1979)
        if alpha == 0
            error('on doit avoir theta different de 0 pour le copule de Franck');
        end
        y(nonzero) = -log((exp(-alpha.*x(nonzero))-1)./(exp(-alpha)-1)); 
        
    case 'joe' %Copule de Joe (Joe 1997)
        if alpha < 1
            error('on doit avoir theta >= 1 pour le copule de Joe');
        end
        
        y(nonzero) = -log(1-(1-x(nonzero)).^alpha);  
        
    
otherwise
    error('family de copule archimédien inconnu: ',family);
end