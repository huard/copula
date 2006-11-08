function d = dilog(z)
% DILOG = di-Logarithm.
%
% d = dilog(z) = Li_2(z) 
%   = -Int From t=0 To t=z    log(1-t) dt/t         for all z.
%   =  Sum From n=1 To n=Inf  z^n/n^2               for |z|<=1.
%
% INPUT  z: real or complex, scalar, vector or matrix.
% OUTPUT d: component-wise dilogarithm of z.

% References:
% [1] Lewin, L. 1958. Dilogarithms and associated functions. Macdonald.
% [2] Wood, D. C. 1992. Technical Report 15-92. University of Kent computing laboratory.
% [3] http://en.wikipedia.org/wiki/Polylog

% Didier Clamond, February 28th, 2006.

% Initialization.
d  = zeros(size(z));     
s  =  ones(size(z));     
      
% For large moduli: Mapping onto the unit circle |z|<=1.
j = find(abs(z)>1);
d(j) = -1.64493406684822643 - 0.5*log(-z(j)).^2; 
s(j) = -s(j);
z(j) = 1./z(j); 

% For large positive real parts: Mapping onto the unit circle with Re(z)<=1/2.
j = find(real(z)>0.5);
d(j) = d(j) + s(j).*( 1.64493406684822643 - log((1-z(j)).^log(z(j))) );
s(j) = -s(j);
z(j) = 1 - z(j);

% Transformation to Debye function and rational approximation.
z = -log(1-z);                                                                                
s = s.*z;                                                                                      
d = d - 0.25*s.*z;                                                                            
z = z.*z;                                                                                     
s = s.*(1+z.*(6.3710458848408100e-2+z.*(1.04089578261587314e-3+z*4.0481119635180974e-6)));    
s = s./(1+z.*(3.5932681070630322e-2+z.*(3.20543530653919745e-4+z*4.0131343133751755e-7)));    
d = d + s;                                                                                    