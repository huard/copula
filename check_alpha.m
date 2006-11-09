function boolean = check_alpha(family, alpha)
%
%   FUNCTION BOOLEAN = CHECK_ALPHA(FAMILY, ALPHA [, WARN])
%   
%   Check ALPHA is a valid parameter for the copula family.
%
%   INPUT
%       FAMILY: One of {'amh' 'arch12' 'arch14' 'clayton' 'frank'
%               'gaussian' 't' 'fgm' 'gumbel'}
%       ALPHA:    Array of copula parameters. 
%   
%   OUTPUT
%       BOOLEAN: Boolean array. True if tau is in the domain, False otherwise.
%       

% D. Huard, Nov. 2006

switch lower(family)
    
    case {'gaussian' 't' 'fgm'}
        boolean = abs(alpha)<=1;
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
        case {'gaussian' 't' 'fgm'}      
            warning('COPULA:BadParameter', 'ALPHA must be in [-1, 1] for the %s copula.\nBad parameters: %s ', family, wrong);
            
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
