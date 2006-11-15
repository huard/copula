function p = posterior(par,data,type)
% P = posterior(par, data, type)
% Compute probability of copula family with given parameter, knowing the data.
% Input:
% par: Copula parameter
% data: NX2 vector (u,v)
% type: Copula family
% G. Evin, 2006

% Compute density
d = copulapdf(data(:,1),data(:,2),type,par);
 
% Remove zeros
nz = (d==0);
d(nz) = ones(size(d(nz))).*eps;
y1 = sum(log(d),1);

% Compute prior on the parameter
nz = (prior_par(type,par) ~= 0);
y2(nz) = log(prior_par(type,par(nz)));

% Prior on Tau
[tau_min,tau_max] = int_tau(type);
y3 = log(tau_max-tau_min);
y = y1 + y2 - y3;
p = exp(y);
