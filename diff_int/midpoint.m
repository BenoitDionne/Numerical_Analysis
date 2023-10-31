%  s = midpoint(funct,a,b,m)
%
%  We use the composite midpoint rule to approximate the integral
%    of the function  funct  between  a and  b .  We use
%    2 * m   subdivisions of the interval.
%  The approximation is stored in  s .
%
function s = midpoint(funct,a,b,m)
  N = 2*m;
  h = (b-a)/N;
  if ( m > 1)
    x = linspace(a,b,N+1);
    xm = x(2:2:N);
    s = sum(funct(xm))*(b-a)/m;
  else
    s = funct(a+h)*(b-a);
  end
end
