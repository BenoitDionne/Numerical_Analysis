A = [ [ 1, 0, 0, sqrt(2)/2, 1, 0, 0, 0 ] ;
      [ 0, 1, 0, sqrt(2)/2, 0, 0, 0, 0 ] ;
      [ 0, 0, 1, 0, 0, 0, -1/2, 0 ] ;
      [ 0. 0, 0, -sqrt(2)/2, 0, -1, 1/2, 0 ] ;
      [ 0, 0, 0, 0, -1, 0, 0, 1 ] ;
      [ 0, 0, 0, 0, 0, 1, 0, 0 ] ;
      [ 0, 0, 0, -sqrt(2)/2, 0, 0, sqrt(3)/2, 0 ] ;
      [ 0, 0, 0, 0, 0, 0, -sqrt(3)/2, -1 ] ];
L = zeros(8,8);
for i=1:7
  L(i+1:8,i) = -A(i+1:8,i);
end
U = zeros(8,8);
for i=1:7
  U(i,i+1:8) = -A(i,i+1:8);
end
D = zeros(8,8);
for i=1:8
  D(i,i) = -A(i,i);
end
b = [0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 10000 ; 0 ];
x = [1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1];
tol = 10^(-3);
limit = 100;

disp 'Jacobi iterative method'
Tj = inv(D)*(L+U);
eig(Tj)
jacobi(A,b,x,tol,limit)

disp 'Gauss-Seifel iterative method'
Tgs = inv(D-L)*U;
eig(Tgs)
gausssiedel(A,b,x,tol,limit)

w = 1.2;
disp(sprintf('Relaxation method with omega = %d',w))
Tr = inv(D- w*L)*((1-w)*D+w*U);
eig(Tr)
relaxation(A,b,x,w,tol,limit)
