function pass = check_tau(family, tau, warn)
%
%   FUNCTION PASS = CHECK_TAU(FAMILY, TAU [, WARN])
%
%   Check TAU is a valid Kendall's rank correlation for the copula family.
%
%   INPUT
%       FAMILY: One of {'amh' 'arch12' 'arch14' 'clayton' 'frank'
%               'gaussian' 't' 'fgm' 'gumbel'}
%       TAU:    Array of Kendall's tau. 
%       WARN:   If 1, display a message whenever some TAUs are not valid.
%               If 0, no message is printed. 
%   
%   OUTPUT
%       PASS: Boolean array. True if tau is in the domain, False otherwise.
%       

% D. Huard, Nov. 2006

%alpha = copulaparam(family, tau)
%pass = check_alpha(alpha)

if nargin <= 2
    warn = 1;
end

base = (tau >= -1) & (tau <=1);

switch lower(family)
    
    case 'amh'
        pass = (tau > -0.181726) & (tau < 1/3);
        
    case {'arch12' 'arch14'}
        pass = tau >= 1/3;
        
    case 'clayton'
        pass = tau > 0;
        
    case 'frank'
        pass = tau ~= 0;
        
    case {'gaussian' 't' 'fgm'}
        pass = (tau >= -1) & (tau <=1); 
           
    case 'gumbel'
        pass = tau >= 0;
        
    otherwise
        error('Copula family ''%s'' not recognized.', family)
end

pass = pass & base;

if any(~pass) && (warn ~= 0)
    switch lower(family)
        case {'gaussian' 't' 'fgm'}      
            fprintf('TAU must be in [-1, 1] for the %s copula.', family);
            
        case 'clayton'
            fprintf('TAU must be nonnegative for Clayton copula.');
            
        case {'gumbel'}
            fprintf('TAU must be greater than or equal to 1 for %s copula.', family);
            
        case 'amh'
            fprintf('TAU must be in [-0.181726, 1/3[ for the Ali-Mikhail-Haq copula.');
            
        case 'frank'
            fprintf('TAU must not be equal to 0 for Frank copula.');
                      
        case {'arch12' 'arch14'}
            fprintf('TAU must be greater than or equal to 1/3 for %s copula.', family);
    end
    fprintf('\nInvalid parameters:')
    disp(tau(~pass))
end

