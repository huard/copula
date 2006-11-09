function u2 = condCDFinv(condCDF,u1,p,par, family)
%
%    CONDCDFINV Inverse conditional distribution function
%    
%    U2 = CONDCDFINV(CONDCDF,U1,P,ALPHA) returns U2 such that
%
%      CONDCDF(U1,U2,ALPHA) = P,
%
%  where CONDCDF is a function handle to a function that computes the
%  conditional cumulative distribution function of U2 given U1, for an
%  archimedean copula with parameter ALPHA.
%
% CONDCDFINV uses a simple binary chop search.  Newton's method or the
% secant method would probably be faster.

lower = zeros(size(p));
upper = ones(size(p));
width = 1;
tol = 1e-12;
while width > tol
    u2 = .5 .* (lower + upper);
    lo = feval(condCDF,u1,u2,par, family) < p;
    lower(lo) = u2(lo);
    upper(~lo) = u2(~lo);
    width = .5 .* width;
end
