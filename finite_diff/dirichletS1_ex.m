% Example:
f = @(x,y) x.*y;
gb = @(x) x.*(2-x);
gt = @(x) x.*(2-x).^2;
gl = @(y) y.*(4-y);
gr = @(y) y.*(4-y).^2;
a = 0;
b = 2;
c = 0;
d = 4;
N = 10;
M = 20;

[X,Y,U] = dirichletS1(f,gb,gt,gl,gr,N,M,a,b,c,d);
[XX,YY] = meshgrid(X,Y);

% We need to transpose the matrix U because meshgrid()
% transposes the coordinates.
surf(XX,YY,U')
xlabel('x')
ylabel('y')
zlabel('u')
