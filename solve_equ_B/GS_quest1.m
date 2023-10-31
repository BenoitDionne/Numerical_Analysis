A = [ [ 3, -1, 1 ] ;
      [ 1, 3, -1 ] ;
      [ 2, 1, -4 ] ];
b = [1 ; 1 ; 0];
x = [0 ; 0 ; 0];
tol = 10^(-5);
limit = 100;

disp 'Gauss-Seifel iterative method'
gausssiedel(A,b,x,tol,limit)
