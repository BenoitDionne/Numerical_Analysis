f = @(x) x.*log(x) + (x.^3)/24 - 5*x.^2;
a = 1;
b = 3;
m = 548;

midpoint(f,a,b,m)
