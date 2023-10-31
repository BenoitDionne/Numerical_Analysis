format long
ff = @(t,y) 100*y + 100*t.^2 -2*t -100;
t0 = 0;
t1 = 1;
y0 = 1;

% Modified Euler Method
disp 'Modified Euler Method'
for N = 10:10:40
    [tt,ww] = modified_euler(ff,t0,t1,N,y0);
    disp(sprintf('N = %5d',N));
    disp(sprintf('y(N) = %0.14g',ww(N+1)));
end
    
% Fourth-Order Runge-Kutta Method
disp 'Fourth-Order Runge-Kutta Method'
for N = 10:10:40
    h = (t1-t0)/N;
    [tt,ww] = rgkt4(ff,h,N,t0,y0);
    disp(sprintf('N = %5d', N));
    disp(sprintf('y(N) = %0.14g', ww(N+1)));
end

% Fourth-Order Adams-Bashforth Method
disp 'Fourth-Order Adams-Bashforth Method'
for N = 10:10:40
    [tt,ww] = adams_bashforth_4s(ff,t0,t1,N,y0);
    disp(sprintf('N = %5d', N));
    disp(sprintf('y(N) = %0.14g', ww(N+1)));
end

% Trapezoidal Method
disp 'Trapezoidal Method'
ffp = @(t,y) 100;
max=10;
tol=0.00001;
n=10;
[tt,rr] = trapez(ff,ffp,t0,y0,t1,n,max,tol)

