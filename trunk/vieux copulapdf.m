function c = copulapdf(family, U, theta)

% function c = copulapdf(family, U, theta)
%
% Function that gives the density of a copula for all U = (u,v) and all parameters
%
% INPUTS:   u and v are Nx2 vectors in [0,1]x[0,1]
%           family is one of the list {'ind', 'gaussian', 'gumbel' 'clayton' 'sim' 'frank' 'gb' 'amh' 'joe'}
%           theta is the 1xM vector of copula parameter.
%
% OUTPUTS:  c, a NxM matrix that contains the results of the copula on the
% points (u,v) for the parameters theta 
%
% Guillaume EVIN
%
%  13 May, 2004.

u = U(:,1)
v = U(:,2)
N = size(U)(1)
M = size(theta)(2)
L = size(theta)(1)
if L > 1 && L != N
   error('Number of parameters must be 1, identical to number of couples in U, or a row vector.')

% Shape vectors into matrices to compute the pdf for all values of (u,v)
% and all parameters without looping.
u = repmat(u, 1, M)
v = repmat(v, 1, M)
theta = repmat(theta, N, 1)

%ld = length(u1);
%lp = length(par);

%s1 = size(u1);
%s2 = size(u2);
%s3 = size(par);

% case of the numerical integration with Matlab, must accept a scalar
% in first argument and a vector in second argument
if (length(u1)>1) & (length(u2) == 1)
    u2 = repmat(u2,length(u1),1);
end

%Les vecteurs doivent tous etre en colonne
if s1(2)>1
    u1 = u1';
end
if s2(2)>1
    u2 = u2';
end
if s3(2)>1
    par = par';
end

if length(par)==1
    par = repmat(par,length(u1),1);
end

if (min(s1)>1 | min(s3)>1)|min(s2)>1
    error('la fonction densitecopulaarch n''admet de matrices comme arguments')
end

if length(par)>1 & (length(par)~=length(u1))
    par = repmat(par,1,max(s1));
    u1 = repmat(u1',max(s3),1);
    u2 = repmat(u2',max(s3),1);
    
    par = reshape(par,prod(size(par)),1);
    u1 = reshape(u1,prod(size(u1)),1);
    u2 = reshape(u2,prod(size(u2)),1);
end

% the function gives 0 out [0,1]x[0,1]
c = zeros(size(u1));
nz = (u1>0 & u1<1) & (u2>0 & u2<1);
switch lower(type)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ellipitical copulas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'gaussian'
        rho = par(nz);
        if (rho < -1 | 1 < rho)
            error('RHO must be a correlation coefficient between -1 and 1');
        else
            v1 = norminv(u1(nz));
            v2 = norminv(u2(nz));
            c(nz) = (1./sqrt(1-rho.^2)).*exp(-(v1.^2+v2.^2-(2.*rho).*v1.*v2)./(2*(1-rho.^2)) + (v1.^2+v2.^2)./2);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Archimedean copulas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case {'ind' 'clayton' 'frank' 'frank2' 'gumbel' 'joe' 'arch12' 'arch14_2' 'arch14' 'fgm' 'amh' 'gb'}

        switch lower(type)
            case {'ind'}
                c = ones(size(u1(nz))); % the density for the independant copula is one

            case 'gumbel'
                % the gumbel copula : C(u,v) = exp(-((-ln u)^par(nz) + (-ln
                % v)^par(nz))^(1/par(nz)))
                C = copula(u1(nz),u2(nz) ,'gumbel',par(nz) );
                logC = (-log(C)).^par(nz);
                v1 = ((-log(u1(nz))).^(par(nz)-1))./u1(nz) ;
                v2 = ((-log(u2(nz))).^(par(nz)-1))./u2(nz) ;
                c(nz) = (v1.*v2.*C).*(logC.^((2-2.*par(nz))./par(nz))-(1-par(nz)).*(logC.^((1-2.*par(nz))./par(nz))));

            case 'clayton'
                % the clayton copula : C(u,v) = u1^(-par(nz))+u2^(-par(nz))-1)^(-1/par(nz))
                C1 = copula(u1(nz),u2(nz),'clayton',par(nz));
                c(nz)= (u1(nz).^(-par(nz)-1)).*(u2(nz).^(-par(nz)-1)).*(par(nz)+1).*C1.^(1+2.*par(nz));

            case 'frank'
                % the frank copula : C(u,v) = (-1/par(nz))*(1 +
                % (exp(-par(nz)*u1)-1)(exp(-par(nz)*u2)-1)/(exp(-par(nz))-1))
                a = exp(-par)-1;
                v1 = exp(-par.*(u1));
                v2 = exp(-par.*(u2));
                nz = (((a + (v1-1).*(v2-1)).^2)~=0)&nz;
                v1 = exp(-par(nz).*(u1(nz)));
                v2 = exp(-par(nz).*(u2(nz)));
                a = exp(-par(nz))-1;
                c(nz) = -par(nz).*v1.*v2.*(a./((a + (v1-1).*(v2-1)).^2));
                
            case 'frank_genest'
                % The frank copula from Genest (1987))
 		% The parameterization is different
                v1 = par(nz).^u1(nz) ;
                v2 = par(nz).^u2(nz) ;
                v3 = par(nz).^(u1(nz) +u2(nz));
                c(nz) = ((par(nz)-1).*log(par(nz)).*v3)./(((par(nz)-1)+(v1-1).*(v2-1)).^2);

            case 'joe'
                % the joe copula : C(u,v) = 1 - [(1-u1)^par(nz)+(1-u2)^par(nz)-
                % (1-u1)^par(nz)(1-u2)^par(nz)]^(1/par(nz))
                C = 1 - copula(u1(nz),u2(nz),'joe',par(nz));
                v1=1-u1(nz);
                v2=1-u2(nz);
                c(nz) = (v1.^(par(nz)-1)).*(v2.^(par(nz)-1)).*par(nz).*(C.^(1-par(nz))) + (par(nz)-1).*((v1.^(par(nz)-1)).*(v2.^(par(nz)-1))...
                    .*(1-v1.^par(nz)).*(1-v2.^par(nz))).*(C.^(1-2*par(nz)));

            case 'arch12'
                % the 12th archimedean copula in Nelsen.
                t1 = -1 + u2(nz)  ;
                t2 = 1 ./ u2(nz)  ;
                t3 = t1 .* t2 ;
                t4 = (-t3) .^ par(nz) ;
                t5 = -1 + u1(nz)  ;
                t6 = 1 ./ u1(nz)  ;
                t7 = t5 .* t6 ;
                t8 = (-t7) .^ par(nz) ;
                t9 = t8 + t4 ;
                t11 = 1 ./ par(nz) ;
                t14 = t9 .^ (-2 .* (-1 + par(nz)) .* t11) ;
                t15 = t4 .* t14 ;
                t16 = 3 .* par(nz) ;
                t17 = (-t7) .^ t16 ;
                t19 = 2 .* par(nz) ;
                t20 = (-t7) .^ t19 ;
                t22 = (-t3) .^ t19 ;
                t25 = t8 .* t14 ;
                t26 = (-t3) .^ t16 ;
                t28 = t4 .* t8 ;
                t29 = t9 .^ t11 ;
                t47 = 1 + t29 ;
                t48 = t47 .^ 2 ;
                c(nz) = (t15 .* t17 + 2 .* t14 .* t20 .* t22 + t25 .* t26 - t28 .* t29 + t28 .* ...
                    par(nz) .* t29 + t15 .* par(nz) .* t17 + 2 .* t22 .* t20 .* t14 .* par(nz) + t25 .* par(nz) .* ...
                    t26) .* t6 ./ t5 .* t2 ./ t1 ./ t48 ./ t47 ./ (t20 + 2 .* t28 + t22) ;

            case 'arch14_2'
                % the 14th archimedean copula in Nelsen.
                v1 = (-1 + u1(nz).^(-1/par(nz))).^par(nz);
                v2 = (-1 + u2(nz).^(-1/par(nz))).^par(nz);
                c(nz) = v1.*v2.*((v1+v2).^(-2+1./par(nz))).*(1+(v1+v2).^(1./par(nz))).^(-2-par(nz)).*(-1+par(nz)+2.*par(nz).*...
                    ((v1+v2).^(1./par(nz))))./(par(nz).*u1(nz).*u2(nz).*(-1+u1(nz).^(1./par(nz))).*(-1+u2(nz).^(1./par(nz))));

            case 'arch14'
                % the 14th archimedean copula of the De Matteis thesis : C(u,v)
                t1 = 1./par(nz);
                t2 = u2(nz).^-t1;
                t3 = t2-1;
                t4 = t3.^par(nz);
                t5 = u1(nz).^-t1;
                t6 = t5-1;
                t7 = t6.^par(nz);
                t8 = t7+t4;
                t9 = t8.^t1;
                t10 = 1+t9;
                t12 = t10.^(-par(nz)-2);
                t13 = t4.*t12;
                t17 = t8.^(-2.*(par(nz)-1).*t1);
                t18 = t17.*par(nz);
                t19 = 3.*par(nz);
                t20 = t6.^t19;
                t23 = 2.*par(nz);
                t24 = t6.^t23;
                t25 = t3.^t23;
                t26 = t24.*t25;
                t27 = t12.*t17;
                t31 = t7.*t12;
                t32 = t3.^t19;
                t35 = t7.*t4;
                t37 = t10.^(-par(nz)-1);
                t38 = t37.*t9;
                t52 = u1(nz).^t1;
                t57 = u2(nz).^t1;
                c(nz) = (t13 .* t18 .* t20 + 2 .* t26 .* t27 .* par(nz) + t31 .* t18 .* t32 - t35 ...
                    .* t38 + t35 .* t38 .* par(nz) + t13 .* t17 .* t20 + 2 .* t26 .* t27 + t31 .* ...
                    t17 .* t32) .* t1 ./ u1(nz)  ./ (-1 + t52) ./ u2(nz)  ./ (-1 + t57) ./ (t24 + 2 .* t35 + t25);

            case 'fgm'
                % the Farlie-Gumbel-Morgenstern Copula
                % = ...
                t3 = par(nz) .* u1(nz) ;
                c(nz) = 1 + par(nz) - (2 .* par(nz)) .* u2(nz)  - 2 .* t3 + 4 .* t3 .* u2(nz) ;

            case 'amh'
                % the Ali-Mak-Hail.... Copula
                % = ...

                if par(nz) < -1 | par(nz) >= 1
                    error('par(nz)ameter must be in [-1, 1) for the AMH copula.');
                else
                    t2 = par(nz) .* u1(nz) ;
                    t3 = t2 .* u2(nz) ;
                    t4 = par(nz).^2;
                    t5 = t4.*u1(nz) ;
                    t7 = par(nz).*u2(nz) ;
                    t10 = -1 + par(nz) - t7 - t2 + t3;
                    t11 = t10.^2;
                    c(nz) = -(1 - 2 .* par(nz) + t3 + t5 .* u2(nz)  + t2 + t7 - t5 - t4 .* u2(nz)  + t4) ./ t11 ./ t10;
                end


            case 'gb'
                % gives the density using the formula with the generator of the archimedean
                % copulas
                g2 = generateurarchderivee2(copula(u1(nz) ,u2(nz) ,type,par(nz)),type,par(nz));
                gu1 = generateurarchderivee(u1(nz) ,type,par(nz));
                gu2 = generateurarchderivee(u2(nz) ,type,par(nz));
                g1 = generateurarchderivee(copula(u1(nz) ,u2(nz) ,type,par(nz)),type,par(nz));
                c(nz) = -g2.*gu1.*gu2./((g1.^3)+(g1.^3==0));

            case 'sort'
            otherwise
                error('Type de copule archimédien inconnu: ',type);
        end
end

