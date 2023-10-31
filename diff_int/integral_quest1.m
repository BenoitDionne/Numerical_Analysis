f = @(x) (x.^2).*log(x);
a = 1;
b = 3;
m = 417;
s0 = 9*log(3) -26/9;

s = midpoint(f,a,b,m)
abs(s-s0)

m = 589;
s = trapezoid(f,a,b,m)
abs(s - s0)

m = 7;
s = simpson(f,a,b,m)
abs(s - s0)
