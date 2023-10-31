p = [ 0 , 0.3 , 1 ];
% f = [ exp(-2*0), exp(-2*0.3), exp(-2*1) ]
% fp = [ -2*exp(-2*0) , -2*exp(-2*1) ]
f = [ 1 , 0.548811636094027 , 0.135335283236613 ];
fp = [ -2 , -0.270670566473225 ];
[L,D,U,b,d] = clampedsplinematrix(f,fp,p)

[z,d] = tridmatrix(L,D,U,b)

x = 0:0.01:1;
[y,coeffs,d] = splinepoly(z,f,p,x);
coeffs
plot(x,y,'b')
grid on
