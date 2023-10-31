format long

f = @(t) [ 0 ; -3*exp(t) ];
A = @(t) [ 0 1 ; 4 0 ];
Ba = [ 1 0 ; 0 0 ];
Bb = [ 0 0 ; 1 0 ];
yc = [ 1 ; exp(1) ];
N = 25;
[t,w] = shooting(f,A,Ba,Bb,yc,N,0,1);

M = [t' w'];
save results.dat M -ascii

% The exact solution
t = linspace(0,1,N+1);
y = [exp(t) ; exp(t) ];
w - y
