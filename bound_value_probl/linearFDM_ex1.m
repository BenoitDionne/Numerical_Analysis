format long

F = @(t,h) [ 0 ; -3*(exp(t)+exp(t+h))/2 ];
Ci = @(t,h) [ -1/h -1/2 ; -2 -1/h];
Cii = @(t,h) [ 1/h -1/2 ; -2  1/h];
Baa = [ 1 0 ];
Bab = [ 0 0 ];
Bbb = [ 1 0 ];
yc = [ 1 ; exp(1) ];

N = 100;
[t,w] = linearFDM(F,Ci,Cii,Baa,Bab,Bbb,yc,N,0,1);

% MM = [t' w'];
% save results.dat MM -ascii -double

% The exact solution
t = linspace(0,1,N+1);
y = exp(t);
w(1,:) - y
