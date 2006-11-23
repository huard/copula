function p = copCondCDF(family, u1,u2,alpha)
%
%   P = COPCONDCDF (FAMILY, U1, U2, ALPHA)
%
%   Conditional cumulative distribution function C(U2|U1, ALPHA)
%
%   INPUT
%	    FAMILY: One of {'GB' (Gumbel & Barnett), Gumbel, Joe, AMH, Arch12,
%         	    Arch14, Clayton, FGM and Frank.}
%	    U1:     Quantile.
%	    U2: 	Quantile.
%	    ALPHA:	Copula parameter.    
%
%   OUTPUT
%	    P:  	Conditional cumulative distribution of U2, given U1.
%

%   G. Evin, 2006

p = zeros(size(u1)); 
nz = ((u1~=0) & (u2~=0));
u1=u1(nz);
u2=u2(nz);

switch lower(family)
    case 'gb'
        p(nz) = (u2+u1.*u2.*(-alpha.*log(u2)./u1)).*exp(-alpha.*log(u1).*log(u2));
    case 'gumbel'
        nlog1 = -log(u1);
        nlog2 = -log(u2);
        C = copulacdf('gumbel', [u1 u2],alpha);
        p(nz) = (nlog1.^(alpha-1)./u1).*((nlog1.^alpha + nlog2.^alpha).^(1/alpha-1)).*C;
    case 'joe'
        mu1 = 1-u1;
        mu2 = 1-u2;
        p(nz) = (mu1.^(alpha-1)).*(1-mu2.^alpha).*((mu1.^alpha+mu2.^alpha-(mu1.^alpha).*(mu2.^alpha)).^(1/alpha-1));
    case 'amh'
        mu1 = 1-u1;
        mu2 = 1-u2;
        p(nz) = ((u2.*(1-alpha.*mu1.*mu2))-u1.*u2.*alpha.*mu2)./((1-alpha.*mu1.*mu2).^2);
    case 'arch12'
        t1 = -1 + u1;
        t2 = 1 ./ u1;
        t4 = (-t1 .* t2) .^ alpha;
        t8 = (-(-1 + u2) ./ u2) .^ alpha;
        t9 = t4 + t8;
        t10 = 1 ./ alpha;
        t11 = t9 .^ t10;
        t13 = (1 + t11) .^ 2;
        t17 = t9 .^ (-(-1 + alpha) .* t10);
        p(nz) = -1 ./ t13 .* t17 .* t4 .* t2 ./ t1;
    case 'arch14'
        t1 = 1 ./ alpha;
        t2 = u1 .^ (-t1);
        t3 = t2 - 1;
        t4 = t3 .^ alpha;
        t5 = u2 .^ (-t1);
        t7 = (t5 - 1) .^ alpha;
        t8 = t4 + t7;
        t9 = t8 .^ t1;
        t11 = -alpha - 1;
        t12 = (1 + t9) .^ t11;
        t13 = alpha - 1;
        t15 = t8 .^ (-t13 * t1);
        t17 = t3 .^ t13;
        t19 = u1 .^ (t11 .* t1);
        p(nz) = t12 .* t15 .* t17 .* t19;
    case 'clayton'
        p(nz) = u1.^(-alpha-1).*(u1.^(-alpha)+u2.^(-alpha)-1).^(-1/alpha-1);
    case 'fgm'
        p(nz) = u2.*(1+alpha.*(1-u2).*(1-2.*u1));
    case 'frank'
        v2 = exp(-alpha.*u2)-1;
        v1 = exp(-alpha.*u1)-1;
        p(nz) = (v2.*(v1+1))./(exp(-alpha)-1+v1.*v2);
end
