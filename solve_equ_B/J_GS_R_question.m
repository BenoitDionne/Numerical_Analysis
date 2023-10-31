A = [ [ 4, 1, -1, 1 ] ;
      [ 1, 5, 1, -1 ] ;
      [ 2, -1, 6, 2 ] ;
      [ -1. 1, -2, 5 ] ];
b = [5 ; 2 ; -14 ; 25];
x = [1 ; 1 ; 1 ; 1];
tol = 10^(-5);
limit = 100;

disp 'Jacobi iterative method'
jacobi(A,b,x,tol,limit)

disp 'Gauss-Seifel iterative method'
gausssiedel(A,b,x,tol,limit)

w = 0.96;
disp(sprintf('Relaxation method with omega = %d', w))
relaxation(A,b,x,w,tol,limit)

% To prove that the relaxation iterative method converges
U = zeros(4,4);
U(1,2:4) = -A(1,2:4);
U(2,3:4) = -A(2,3:4);
U(3,4) = -A(3,4);
L = zeros(4,4);
L(2:4,1) = -A(2:4,1);
L(3:4,2) = -A(3:4,2);
L(4,3) = -A(4,3);
D = zeros(4,4);
D(1,1) = A(1,1);
D(2,2) = A(2,2);
D(3,3) = A(3,3);
D(4,4) = A(4,4);
T = inv(D-w*L)*((1-w)*D+w*U)
% c = w*inv(D-w*L)*b;
eig(T)
