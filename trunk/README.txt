COPULA BCS
Hosted at code.google.com/p/copula/

This set of Matlab (TM) files provides the functions needed to estimate the 
best copula, given a set of copulas and fractiles. 

The method is described in the paper
Huard, D., Évin, G. and Favre, A-C. Bayesian Copula Selection, 
Journal of Computational Statistics and Data Analysis, 2005, 51, 809-822.

Currently, only a subset of bivariate copulas are supported, namely
Clayton, Gumbel, Frank, Gaussian, AMH, FGM, Arch12, Arch14, Ind.


The main function is bcs.m, which computes the probability for a collection of copula families. 
See tests/example.m to see how it works. 


Functions description

--------------------


bcs: Return the weight of each copula family, given the data.

check_alpha: Return a boolean indicating the validity of the copula parameter.

check_tau: Return a boolean indicating the validity of Kendall's tau for a given copula.

constrain_tau: Constrain bounds in tau to the domain convered by the copula.  

copulacdf: Return the cumulative distribution function of the copula.

copula_like: Return the likelihood of data set.

copulaparam: Find copula parameter given Kendall's tau.

copulapdf: Return the probability distribution function of the copula.

copulastat: Kendall's rank correlation for a copula.

tauboundaries: Return the domain spanned by Kendall's tau.

taujacobian: Return the derivative of tau(alpha) with respect to alpha.

TODO: Find a solution to finite precision problems, ie NaN and Infs creeping out from copulapdf. 

The code is written mostly by G. Evin and D. Huard. It extends functions originaly written by P. Perkins.

To suggest improvements or to signal bugs, use the Issue tracker at 
http://code.google.com/p/copula/issues/list


David Huard 
November 22, 2006


