function c = copulapdf(family, U, alpha)
%
%   FUNCTION C = COPULAPDF(FAMILY, U, ALPHA)
%
%   Return c(u,v|ALPHA), the copula density. 
%
%   INPUTS
%       FAMILY: one of {'ind', 'gaussian', 'gumbel' 'clayton' 'frank' 'amh' 'joe' 'fgm' 'arch12' 'arch14'}       
%       U: Nx2 vector (u,v) in [0,1]^2.
%       ALPHA: 1xM vector of copula parameters
%           or 
%       U: 1x2 vector (u,v) in [0,1]^2.
%       ALPHA: NxM vector of copula parameters. 
%
%   OUTPUTS
%       C: NxM matrix of C(u,v|ALPHA) 

%   Guillaume EVIN, 13 May, 2004.
%   D. Huard, Nov. 2006
%  

% Check alpha is in the domain covered by the family.
pass = check_alpha(family, alpha);
if ~all(pass)
    error('Some parameters are not valid.\n%f', alpha(~pass))
end

% Check u,v are in [0,1]^2
if any( (U < 0) | (U > 1) )
    error('Some quantiles are outside the unit hypercube.')
end

% Shape checking
[NU, MU] = size(U);
[NA, MA] = size(alpha);

if MU ~= 2
    error('Bad shape. U is not Nx2, but rather %s.', mat2str(size(U)))
end

% Reshape ALPHA
if NA == 1 
    alpha = repmat(alpha, NU, 1);
elseif NA ~= NU && NU ~= 1
    error('Number of parameters must be 1, identical to number of couples in U, or a row vector.')
end

% Reshape u,v
u = U(:,1);
v = U(:,2);
if NU == 1
    u = repmat(u, NA, MA);
    v = repmat(v, NA, MA);
else
    u = repmat(u, 1, MA);
    v = repmat(v, 1, MA);
end

switch lower(family)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ellipitical copulas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 'gaussian'
        v1 = norminv(u);
        v2 = norminv(v);
        c = (1./sqrt(1-alpha.^2)).*exp(-(v1.^2+v2.^2-(2.*alpha).*v1.*v2)./(2*(1-alpha.^2)) + (v1.^2+v2.^2)./2);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Archimedean copulas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'ind'
        c = ones(size(u)); % the density for the independant copula is one
        
    case 'gumbel'
        % the gumbel copula : C(u,v) = exp(-((-ln u)^alpha + (-ln
        % v)^alpha)^(1/alpha))
        C = copulacdf('gumbel',[u(:,1),v(:,1)], alpha(1,:));
        logC = (-log(C)).^alpha;
        v1 = ((-log(u)).^(alpha-1))./u ;
        v2 = ((-log(v)).^(alpha-1))./v ;
        c = (v1.*v2.*C).*(logC.^((2-2.*alpha)./alpha)-(1-alpha).*(logC.^((1-2.*alpha)./alpha)));
        
    case 'clayton'
        % the clayton copula : C(u,v) = u^(-alpha)+v^(-alpha)-1)^(-1/alpha)
        C1 = copulacdf('clayton',[u(:,1),v(:,1)], alpha(1,:));
        c= (u.^(-alpha-1)).*(v.^(-alpha-1)).*(alpha+1).*C1.^(1+2.*alpha);
        
    case 'frank'
        % the frank copula : C(u,v) = (-1/alpha)*(1 +
        % (exp(-alpha*u)-1)(exp(-alpha*v)-1)/(exp(-alpha)-1))
        a = exp(-alpha)-1;
        v1 = exp(-alpha.*(u));
        v2 = exp(-alpha.*(v));
        %        nz = (((a + (v1-1).*(v2-1)).^2)~=0);
        v1 = exp(-alpha.*(u));
        v2 = exp(-alpha.*(v));
        a = exp(-alpha)-1;
        c = -alpha.*v1.*v2.*(a./((a + (v1-1).*(v2-1)).^2));
        
    case 'frank_genest'
        % The frank copula from Genest (1987))
        % The alphaameterization is different
        v1 = alpha.^u ;
        v2 = alpha.^v ;
        v3 = alpha.^(u +v);
        c = ((alpha-1).*log(alpha).*v3)./(((alpha-1)+(v1-1).*(v2-1)).^2);
        
    case 'joe'
        % the joe copula : C(u,v) = 1 - [(1-u)^alpha+(1-v)^alpha-
        % (1-u)^alpha(1-v)^alpha]^(1/alpha)
        C = 1 - copulacdf('joe',[u(:,1),v(:,1)], alpha(1,:));
        v1=1-u;
        v2=1-v;
        c = (v1.^(alpha-1)).*(v2.^(alpha-1)).*alpha.*(C.^(1-alpha)) + (alpha-1).*((v1.^(alpha-1)).*(v2.^(alpha-1))...
            .*(1-v1.^alpha).*(1-v2.^alpha)).*(C.^(1-2*alpha));
        
    case 'arch12'
        % the 12th archimedean copula in Nelsen.
        t1 = -1 + v  ;
        t2 = 1 ./ v  ;
        t3 = t1 .* t2 ;
        t4 = (-t3) .^ alpha ;
        t5 = -1 + u  ;
        t6 = 1 ./ u  ;
        t7 = t5 .* t6 ;
        t8 = (-t7) .^ alpha ;
        t9 = t8 + t4 ;
        t11 = 1 ./ alpha ;
        t14 = t9 .^ (-2 .* (-1 + alpha) .* t11) ;
        t15 = t4 .* t14 ;
        t16 = 3 .* alpha ;
        t17 = (-t7) .^ t16 ;
        t19 = 2 .* alpha ;
        t20 = (-t7) .^ t19 ;
        t22 = (-t3) .^ t19 ;
        t25 = t8 .* t14 ;
        t26 = (-t3) .^ t16 ;
        t28 = t4 .* t8 ;
        t29 = t9 .^ t11 ;
        t47 = 1 + t29 ;
        t48 = t47 .^ 2 ;
        c = (t15 .* t17 + 2 .* t14 .* t20 .* t22 + t25 .* t26 - t28 .* t29 + t28 .* ...
            alpha .* t29 + t15 .* alpha .* t17 + 2 .* t22 .* t20 .* t14 .* alpha + t25 .* alpha .* ...
        t26) .* t6 ./ t5 .* t2 ./ t1 ./ t48 ./ t47 ./ (t20 + 2 .* t28 + t22) ;
        
    case 'arch14'
        % the 14th archimedean copula in Nelsen.
        v1 = (-1 + u.^(-1./alpha)).^alpha;
        v2 = (-1 + v.^(-1./alpha)).^alpha;
        c = v1.*v2.*((v1+v2).^(-2+1./alpha)).*(1+(v1+v2).^(1./alpha)).^(-2-alpha).*(-1+alpha+2.*alpha.*...
            ((v1+v2).^(1./alpha)))./(alpha.*u.*v.*(-1+u.^(1./alpha)).*(-1+v.^(1./alpha)));
        
    case 'fgm'
        % the Farlie-Gumbel-Morgenstern Copula
        % = ...
        t3 = alpha .* u ;
        c = 1 + alpha - (2 .* alpha) .* v  - 2 .* t3 + 4 .* t3 .* v ;
        
    case 'amh'
        % the Ali-Mikhail-Haq Copula
        % = ...
        t2 = alpha .* u ;
        t3 = t2 .* v ;
        t4 = alpha.^2;
        t5 = t4.*u ;
        t7 = alpha.*v ;
        t10 = -1 + alpha - t7 - t2 + t3;
        t11 = t10.^2;
        c = -(1 - 2 .* alpha + t3 + t5 .* v  + t2 + t7 - t5 - t4 .* v  + t4) ./ t11 ./ t10;
        
        
        
        %         case 'gb'
        %             % gives the density using the formula with the generator of the archimedean
        %             % copulas
        %             g2 = generateurarchderivee2(copula(u ,v ,type,alpha),type,alpha);
        %             gu = generateurarchderivee(u ,type,alpha);
        %             gv = generateurarchderivee(v ,type,alpha);
        %             g1 = generateurarchderivee(copula(u ,v ,type,alpha),type,alpha);
        %             c = -g2.*gu.*gv./((g1.^3)+(g1.^3==0));
        
    otherwise
        error('Copula family ''%s'' not recognized.', family)
end
