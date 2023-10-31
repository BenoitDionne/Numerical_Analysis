err = sqrt((exp(1)*10^4)/4)
m = ceil(err)
f = @(x) exp(x.^2);
midpoint(f,0,1,m)
