format long

f = @(t) [ 0 ; 0 ; 0 ; -1 + 200*t.^2 ];
A = @(t) [ 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1; -400 0 401 0];
Ba = [ 1 0 0 0 ; 0 1 0 0 ; 0 0 0 0 ; 0 0 0 0 ];
Bb = [ 0 0 0 0 ; 0 0 0 0 ; 1 0 0 0 ; 0 1 0 0 ];
yc = [ 1 ; 1 ; 3/2 + sinh(1) ; 1 + cosh(1) ];

N = 25;
[t,w] = shooting(f,A,Ba,Bb,yc,N,0,1);

% M = [t' w(1,:)'];
% save results.dat M -ascii -double

% The exact solution
t = linspace(0,1,N+1);
y = 1 + t.^2/2 + sinh(t);
w(1,:) - y

% plot(t,y,'k')
% hold on
% plot(t,w(1,:),'b')
% xlabel('t')
% ylabel('y')
% grid on

% funct1 = @(t,y) A(t)*y + f(t);
% N = 20;
% h = 1/N;
% t0 = 0;
% [t,w] = rgkt4(funct1,h,N,t0,yc)

