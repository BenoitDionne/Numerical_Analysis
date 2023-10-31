%  function euler_ex(funct, sol, t0, tf, N, x0)
%    funct is the right hand side  in x' = f(t,x)
%    sol is the solution for the initial condition x(t_0) = x_0
%    t0 is t_0
%    tf is t_f
%    N is the number of intervals
%    x0 is x_0
%
% First example:
% ff = @(t,x) 0.2*t*x;
% ss = @(t) exp(0.1*t.^2 -0.1);
% t0=1, t1=1.5, N=5, x0 = 1;
% result = euler_ex(ff, ss, t0, t1, N, x0)
% save results.dat result -ascii
%
% Second example:
% ff = @(t,x) 2*t*x;
% ss = @(t) exp(t.^2 - 1);
% t0=1, t1=1.5, N=5, x0 = 1;
% result = euler_ex(ff, ss, t0, t1, N, x0)
% save results.dat result -ascii
%
% Third example:
% ff = @(t,x) 1 - 2*x;
% ss = @(t) (1+exp(-2*t))/2;
% t0=0; t1=4; N=128; x0 = 1;
% result = euler_ex(ff, ss, t0, t1, N, x0)
% save results.dat result -ascii
%
function result = euler_ex(funct, sol, t0, t1, N, x0)
  result = zeros(N+1,6);
  h = (t1-t0)/N;
  result(1,1) = 0;
  result(1,2) = t0;
  result(1,3) = x0; 
  result(1,4) = sol(t0);
  result(1,5) = 0;
  result(1,6) = 0;
  for n=1:1:N
      ni = n+1;
      result(ni,1) = n;
      result(ni,2) = t0 + n*h;
      result(ni,3) = result(n,3) + h.*funct(result(n,2),result(n,3));
      result(ni,4) = sol(result(ni,2));
      result(ni,5) = result(ni,3) - result(ni,4);
      result(ni,6) = abs(result(ni,5)/result(ni,4));
  end
end
