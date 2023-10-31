%  function adams_bashforth_2s(funct, sol, t0, tf, N, x0)
%    funct is the right hand side  in x' = f(t,x)
%    sol is the solution for the initial condition x(t_0) = x_0
%    t0 is t_0
%    tf is t_f
%    N is the number of intervals
%    x0 is x_0
%
% Example:
% ff = @(t,x) 1 - 2*x;
% ss = @(t) (1+exp(-2*t))/2;
% t0=0; t1=4; N=128; x0 = 1;
% result = adams_bashforth_2s(ff, ss, t0, t1, N, x0)
%
function result = adams_bashforth_2s(funct, sol, t0, t1, N, x0)
  result = zeros(N+1,6);
  h = (t1-t0)/N;
  h2 = 2*h;
  result(1,1) = 0;
  result(1,2) = t0;
  result(1,3) = x0; 
  result(1,4) = sol(t0);
  result(1,5) = 0;
  result(1,6) = 0;  
  result(2,1) = 1;
  result(2,2) = t0+h;
  result(2,3) = sol(result(2,2)); 
  result(2,4) = result(2,3);
  result(2,5) = 0;
  result(2,6) = 0;
  for n=2:1:N
      ni = n+1;
      result(ni,1) = n;
      result(ni,2) = t0 + n*h;
      result(ni,3) = result(n-1,3) + h2*funct(result(n,2),result(n,3));
      result(ni,4) = sol(result(ni,2));
      result(ni,5) = result(ni,3) - result(ni,4);
      result(ni,6) = abs(result(ni,5)/result(ni,4));
  end
end
