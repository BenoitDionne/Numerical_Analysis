format long
ff = @(t,y) 1 + (t-y).^2;
ss = @(t) t - 1/(t-1);
t0=2; t1=3; N=10; y0 = 1;
result = runge_kutta_ex(ff, ss, t0, t1, N, y0)
save results.dat result -ascii
