function boolean = bool_par(par,type)

switch lower(type)
    case {'gaussian' 't' 'fgm'}
        boolean = abs(par)<=1;
    case 'clayton'
        boolean = (par >= 0);
    case 'frank'
        boolean = (par ~= 0);
    case 'gb'
        boolean = (par >= 0)&(par<=1);
    case {'gumbel' 'joe' 'arch12' 'arch14'}
        boolean = (par >= 1);
    case 'amh'
        boolean = (par >= -1)&(par<1);
    otherwise
        error('unknown copula family.');
end