format long

f = @(t) [ 0 ; -3*exp(t) ];
A = @(t) [ 0 1 ; 4 0 ];
Baa = [ 1 0 ];
Bab = [ 0 0 ];
Bbb = [ 1 0 ];
yc = [ 1 ; exp(1) ];

N = 10;
M = 10; 
[t,w] = par_shooting(f,A,Baa,Bab,Bbb,yc,M,N,0,1);

MM = [t' w'];
save results.dat MM -ascii

% The exact solution
t = linspace(0,1,N*M+1);
y = [exp(t) ; exp(t) ];
w - y
