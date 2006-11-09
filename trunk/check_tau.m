function boolean = check_tau(family, tau)
%
%   FUNCTION BOOLEAN = CHECK_TAU(FAMILY, TAU)
%
%   Check TAU is a valid Kendall's rank correlation for the copula family.
%
%   INPUT
%       FAMILY: One of {'amh' 'arch12' 'arch14' 'clayton' 'frank'
%               'gaussian' 't' 'fgm' 'gumbel'}
%       TAU:    Array of Kendall's tau. 
%          
%   OUTPUT
%       BOOLEAN: Boolean array. True if tau is in the domain, False otherwise.
%       

%   G. Evin, D. Huard, Nov. 2006

%alpha = copulaparam(family, tau)
%pass = check_alpha(alpha)

base = (tau >= -1) & (tau <=1);

switch lower(family)
    
    case {'gaussian' 't' 'plackett'}
        boolean = ones(size(tau));
    case 'clayton'
        boolean = (tau > 0);
    case 'frank'
        boolean = (tau ~= 0);
    case 'fgm'
        boolean = (abs(tau) <= 2/9);
    case 'gb'
        boolean = (tau <= 0);
    case 'gumbel'
        boolean = (tau >= 0);
    case 'amh'
        boolean = (tau >= -0.181726)*(tau <= 1/3);
    case 'joe'
        boolean = (tau > 0);
    case {'arch12' 'arch14'}
        boolean = (tau >= 1/3);
    otherwise
        error('Copula family ''%s'' not recognized.', family)
end

boolean = boolean & base;

if any(~boolean) && (warn ~= 0)
    wrong = mat2str(tau(~boolean));
    switch lower(family)
        case {'gaussian' 't' 'plackett'}      
            warning('COPULA:BadParameter', 'TAU must be in [-1, 1] for the %s copula.\nBad parameters: %s', family, wrong);
            
        case 'clayton'
            warning('COPULA:BadParameter', 'TAU must be in ]0,1] for Clayton copula.\nBad parameters: %s',wrong);
            
        case 'frank'
            warning('COPULA:BadParameter', 'TAU must be in [-1,1]\{0} for Frank copula.\nBad parameters: %s',wrong);

        case {'fgm'}      
            warning('COPULA:BadParameter', 'TAU must be in [-2/9, 2/9] for FGM copula.\nBad parameters: %s',wrong);

        case {'gb'}      
            warning('COPULA:BadParameter', 'TAU must be in [-1,0] for Gumbel-Barnett copula.\nBad parameters: %s',wrong);

        case {'gumbel'}
            warning('COPULA:BadParameter', 'TAU must be in [0,1] for Gumbel copula.\nBad parameters: %s',wrong);
            
        case 'amh'
            warning('COPULA:BadParameter', 'TAU must be in [-0.181726, 1/3] for Ali-Mikhail-Haq copula.\nBad parameters: %s',wrong);
                    case 'joe'
            warning('COPULA:BadParameter', 'TAU must be in ]0, 1] for Ali-Mikhail-Haq copula.\nBad parameters: %s',wrong);

        case {'arch12' 'arch14'}
            warning('COPULA:BadParameter', 'TAU must be in [1/3,1] for %s copula.\nBad parameters: %s', family, wrong);
    end
end

