%  s = trapezoid(funct,a,b,N)
%
%  We use composite trapezoidal rule to approximate the integral
%  of a function.
%
%  The following must be given:
%    The function given in  "funct",
%    The endpoints  a  and  b  of the interval  [a,b],
%    The number  N  of subintervals to be used.
%
%  The program gives an approximation  s  of the integral.
%
function s = trapezoid(funct,a,b,N)
  h = (b-a)/N;
  sp = (funct(a) + funct(b))/2;
  if ( N > 1)
    x = linspace(a,b,N+1);
    x2 = x(2:N);
    s = (sp + sum(funct(x2)))*h;
  else
    s = sp*h;
  end
end
