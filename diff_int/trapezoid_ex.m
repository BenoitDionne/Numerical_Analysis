err = sqrt((exp(1)*10^4)/2)
n = ceil(err)
f = @(x) exp(x.^2);
trapezoid(f,0,1,n)

