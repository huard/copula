function y=frankdg(x)
% Computes the derivative of g(theta) where Kendall's tau=g(theta) and theta is the parameter of the Frank Copula.;
t2 = exp(x); 
t5 = 1-t2;
t6 = log(t5); 
t7 = x.*t6; 
s11 = maple('a:=polylog',2,t2);
s12 = maple('evalf(a)');
t11 = str2num(s12);
t15 = pi^2; 
t17 = x^2; 
t27 = real(-(-3*x+3.*x.*t2+6*t7-6*t7.*t2+6*t11-6*t2.*t11-t15+t2.*t15+3*t17.*t2)./t17/x./t5); 
y = 4/3*t27; 