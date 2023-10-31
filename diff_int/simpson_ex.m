err = ( (exp(1)*10^4*76)/(180*2^4) )^(1/4)
m = ceil(err)
f = @(x) exp(x.^2);
simpson(f,0,1,m)
