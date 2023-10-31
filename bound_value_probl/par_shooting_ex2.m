format long

f = @(t) [ 0 ; 0 ; 0 ; -1 + 200*t.^2 ];
A = @(t) [ 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1; -400 0 401 0];
Baa = [ 1 0 0 0 ; 0 1 0 0 ];
Bab = [ 0 0 0 0 ; 0 0 0 0 ];
Bbb = [ 1 0 0 0 ; 0 1 0 0 ];
yc = [ 1 ; 1 ; 3/2 + sinh(1) ; 1 + cosh(1) ];

M = 10;
N = 10;
[t,w] = par_shooting(f,A,Baa,Bab,Bbb,yc,M,N,0,1);

MM = [t' w(1,:)'];
save results.dat MM -ascii -double

% The exact solution
t = linspace(0,1,M*N+1);
y = 1 + t.^2/2 + sinh(t);
w(1,:) - y

% plot(t,y,'k')
% hold on
% plot(t,w(1,:),'b')
% xlabel('t')
% ylabel('y')
% grid on
