% Since the rasulting matrix A in singular, the trapezoidal method
% cannot be used for this problem.  So, this program fails.

format long

F = @(t,h) [ 0 ; 0 ; 0 ; -1 + 200*t.^2 ];
Ci = @(t,h) [ -1/h -1/2 0 0 ; 0 -1/h -1/2 0 ; 0 0 -1/h -1/2; 200 0 -401/2 -1/h];
Cii = @(t,h) [ 1/h -1/2 0 0 ; 0 1/h -1/2 0 ; 0 0 1/h -1/2; 200 0 -401/2 1/h];
Baa = [ 1 0 0 0 ; 0 1 0 0 ];
Bab = [ 0 0 0 0 ; 0 0 0 0 ];
Bbb = [ 1 0 0 0 ; 0 1 0 0 ];
yc = [ 1 ; 1 ; 3/2 + sinh(1) ; 1 + cosh(1) ];

N = 100;
[t,w] = linearFDM(F,Ci,Cii,Baa,Bab,Bbb,yc,N,0,1);

% MM = [t' w(1,:)'];
% save results.dat MM -ascii -double

% The exact solution
% t = linspace(0,1,N+1);
% y = 1 + t.^2/2 + sinh(t);
% w(1,:) - y
