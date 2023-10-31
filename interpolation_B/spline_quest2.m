p = [ 0 , 1 , 3 , 4 , 5 , 5.5 ];
f = [ 2 , 3 , 3 , 2 , 2 , 1 ];
[L,D,U,b] = naturalsplinematrix(f,p)

z = tridmatrix(L,D,U,b);
z = [0,z,0]

x = 0:0.01:5.5;
[y,coeffs] = splinepoly(z,f,p,x);
coeffs
plot(x,y,'b')
grid on
xlabel('x')
ylabel('y')
