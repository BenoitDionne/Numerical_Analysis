%  newton_chaos_quest(a,xl,xu,nmeshes,tol,N)
%
%  For each initial value x0 in the set linspace(xl,xu,nsteps),
%  we use Newton's Method to approximate a real root of the
%  polynomial defined by the vector  a.
%
%  Input:
%    The coefficients a_0 (a(1) in the code below), a_1
%      (a(2) in the code below), ... , a_n (a(max) in the code below),
%      ... of the polynomial
%      p(x) = a_0 + a_1 x + a_2 x^2 + ... + a_n x^n
%      (even the zero coefficients must be given.)
%    The lower bound xl of the interval
%    The upper bound xu of the interval
%    The number of equaly spaced points chosen in the interval [xl,xu]
%    The maximal tolerance  tol
%    The maximal number N of iterations
%  Output:
%    The graph of the (approxiamtions of the) roots of the polynomial
%    as a function of  x0  in Newton's Method.
%
function newton_chaos_quest(a,xl,xu,nmeshes,tol,N)
  max = length(a);
  x0 = linspace(xl,xu,nmeshes);
  roots = [];

  for n=1:nmeshes
    x = x0(n);
    for k=1:N
      y = a(max);
      z = a(max);
      for i=max-1:-1:2
        y = a(i) + x*y;
        z = y + x*z;
      end
      y = a(1) + x*y;
      % y = p(x)  and  z = p'(x) .

      ratio = y/z;
      x = x - ratio;
      if (abs(ratio) < tol)
        roots(n) = x;
        break;
      end
    end
    if (k == N)
      roots(n) = NaN;
    end
  end

  cla
  plot(x0,roots);
  hold on
  xlabel('x_0');
  ylabel('root');
end

