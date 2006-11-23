function boolean = check_alpha(family, alpha)
%
%   FUNCTION BOOLEAN = CHECK_ALPHA(FAMILY, ALPHA])
%   
%   Check ALPHA is a valid parameter for the copula family.
%
%   INPUT
%       FAMILY: One of {'AMH' 'Arch12' 'Arch14' 'Clayton' 'FGM' 'Frank'
%               'Gaussian' 'GB' 'Gumbel' 'Joe' 't'}
%       ALPHA:  Array of copula parameters. 
%               For the t copula, alpha must be [rho, nu].
%   
%   OUTPUT
%       BOOLEAN: Boolean array. True if tau is in the domain, False otherwise.
%       

% G. Evin, D. Huard, Nov. 2006

switch lower(family)
    
    case {'gaussian' 'fgm'}
        boolean = abs(alpha)<=1;
    case 't'
        rho = alpha(1);
        nu = alpha(2);
        boolean = abs(rho)<=1 & nu >= 0 & mod(nu,1)==0;
    case 'clayton'
        boolean = (alpha >= 0);
    case 'frank'
        boolean = (alpha ~= 0);
    case 'gb'
        boolean = (alpha >= 0)&(alpha<=1);
    case {'gumbel' 'joe' 'arch12' 'arch14'}
        boolean = (alpha >= 1);
    case 'amh'
        boolean = (alpha >= -1)&(alpha<1);
    case 'ind'
        boolean = ones(size(alpha));
    otherwise
        error('Copula family ''%s'' not recognized', family)
end

if any(~boolean) 
    wrong = num2str(alpha(~boolean));
    switch lower(family)
        case {'gaussian' 'fgm'}      
            warning('COPULA:BadParameter', 'ALPHA must be in [-1, 1] for the %s copula.\nBad parameters: %s ', family, wrong);
            
        case 't'
            warning('COPULA:BadParameter', 'ALPHA must be in [-1, 1] and DF a positive integer for the t copula.\nBad parameters: %s ', wrong);
        case 'clayton'
            warning('COPULA:BadParameter', 'ALPHA must be nonnegative for the Clayton copula.\nBad parameters: %s', wrong);
            
        case {'gumbel' 'joe'}
            warning('COPULA:BadParameter', 'ALPHA must be greater than or equal to 1 for the %s copula.\nBad parameters: %s', family, wrong);
            
        case 'amh'
            warning('COPULA:BadParameter', 'ALPHA must be in [-1,1[ for the Ali-Mikhail-Haq copula.\nBad parameters: %s', wrong);
            
        case 'frank'
            warning('COPULA:BadParameter', 'ALPHA must not be equal to 0 for Frank copula.\nBad parameters: %s', wrong);
            
        case 'gb'
            warning('COPULA:BadParameter', 'ALPHA must be in [0,1] for the Gumbel-Barnett copula.\nBad parameters: %s', wrong);
            
        case {'arch12' 'arch14'}
            warning('COPULA:BadParameter', 'ALPHA must be greater than or equal to 1 for the %s copula.\nBad parameters: %s', family, wrong);
    end
end
