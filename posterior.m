function p = posterior(par,data,type)
	% P = posterior(par, data, type)
	% Compute probability of copula family with given parameter, knowing the data.
	% Input:
	% par: Copula parameter
	% data: NX2 vector (u,v)
	% type: Copula family

	% Compute density
        d = densitecopulaarch(data(:,1),data(:,2),type,par);

	% Remove zeros
        nz = (d==0);
        d(nz) = ones(size(d(nz))).*eps;

        matdnn=reshape(d,length(par),length(data));
        [tau_min,tau_max] = int_tau(type);
        y1 = sum(log(matdnn),2)';
	
	% Compute prior on the parameter
        nzprior = (prior_par(type,par) ~= 0);
        y2(nzprior) = log(prior_par(type,par(nzprior)));
 	
	% Prior on Tau
        y3 = log(tau_max-tau_min);

        y = y1 + y2 - y3;
	p = exp(y);
