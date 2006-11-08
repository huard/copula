function pass = check_alpha(family, alpha)

%   FUNCTION PASS = CHECK_ALPHA(FAMILY, ALPHA)
%   
%   Check the ALPHA is a valid parameter for the copula family.
%   Print a message if some parameters are not valid.
%

%   D. Huard, Nov. 2006

switch lower(family)
    
    case {'gaussian' 't' 'fgm'}
        pass = (alpha >= -1) & (alpha <=1); 
        
    case 'clayton'
        pass = alpha > 0;
        
    case {'gumbel' 'joe'}
        pass = alpha >= 1;
        
    case 'amh'
        pass = alpha >=-1 & alpha < 1; 
        
    case 'frank'
        pass = alpha ~= 0
        
    case 'gb'
        pass = alpha >=0 & alpha <=1;
        
    case {'arch12' 'arch14'}
        pass = alpha >= 1;
    
    otherwise
        error('Copula family ''%s'' not recognized', family)
end

if any(~pass)
    switch lower(family)
        case {'gaussian' 't'}      
            fprintf('ALPHA must be in [-1, 1] for the %s copula.', family);
            
        case 'clayton'
            fprintf('ALPHA must be nonnegative for the Clayton copula.');
            
        case {'gumbel' 'joe'}
            fprintf('ALPHA must be greater than or equal to 1 for the %s copula.', family);
            
        case 'amh'
            fprintf('ALPHA must be in [-1,1[ for the Ali-Mikhail-Haq copula.');
            
        case 'frank'
            fprintf('ALPHA must not be equal to 0 for Frank copula.');
            
        case 'gb'
            fprintf('ALPHA must be in [0,1] for the Gumbel-Barnett copula.');
            
        case {'arch12' 'arch14'}
            fprintf('ALPHA must be greater than or equal to 1 for the %s copula.', family);
    end
    fprintf('\nInvalid parameters:')
    disp(alpha(~pass))
end