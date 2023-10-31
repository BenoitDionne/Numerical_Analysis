%  s = simpson(funct,a,b,m)
%
%  We use composite Simpson's rule to approximate the integral
%  of a function.
%
%  The following must be given:
%    The function given in  "funct",
%    The endpoints  a  and  b  of the interval  [a,b],
%    The number  m  which is half the number of subintervals to be used.
%
%  The program gives an approximation  s  of the integral.
%
function s = simpson(funct,a,b,m)
  N = 2*m;
  h = (b-a)/N;
  if ( m > 1)
    x = linspace(a,b,N+1);
    x4 = x(2:2:N);
    x2 = x(3:2:N-1);
    s = h*(funct(a) + funct(b) + 2*sum(funct(x2)) + 4*sum(funct(x4)))/3;
  else
    s = h*(funct(a) + funct(b) + 4*funct(a+h))/3;
  end
end
