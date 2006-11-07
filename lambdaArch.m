function u = lambdaarch(t,alpha,type)
%function u = lambdaarch(t,alpha,type)
%
% Function used in the paper of Genest [1993], in order to compare
% differents copulas, lambda(t) = g(t)/g'(t) where g is the generator
%
% INPUTS:   t is a vector that lies in [0,1]
%           type is one of the list {'gumbel' 'clayton' 'frank' 'gb' 'amh' 'joe'}
%           alpha is the parameter of the copula
%
% OUTPUTS:	u, a vectors that contains the results of lambda on the
% points t
%
% Guillaume EVIN
%
%  13 May, 2004.

% avoid a division by zero
den = generateurarchderivee(t,type,alpha);
nonnul = (den~=0);
nul= (den == 0);

u(nonnul) = generateurarch(t(nonnul),type,alpha)./generateurarchderivee(t(nonnul),type,alpha);
u(nul) = sign(generateurarch(t(nul),type,alpha)).*realmax;