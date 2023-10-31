% Test
f = @(x,y) 0;
gb = @(x) 0;
gt = @(x) x.^2/2;
gl = @(y) sin(pi*y);
gr = @(y) exp(pi)*sin(pi*y) +y.^2/2;
a = 0;
b = 1;
c = 0;
d = 1;
N = 20;
M = 20;

[X,Y,U] = dirichletS1(f,gb,gt,gl,gr,N,M,a,b,c,d);

[XX,YY] = meshgrid(X,Y);
surf(XX,YY,U')
xlabel('x')
ylabel('y')
zlabel('u')
