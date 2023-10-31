p = [ 0 , 1 , 2,  3 , 5 ]
% f = [ cos(0) , cos(1), cos(2), cos(3), cos(5) ]
f = [ 1 , 0.540302305868140 , -0.416146836547142, -0.989992496600445 , ...
      0.283662185463226 ];

[L,D,U,b,d] = naturalsplinematrix(f,fp,p)

[z,d] = tridmatrix(L,D,U,b);

z = [0, z(1:3), 0]

x = 0:0.01:5;
[y,coeffs,d] = splinepoly(z,f,p,x);
coeffs

plot(x,y,'b')
hold on
grid on
z = cos(x);
plot(x,z,'k')
