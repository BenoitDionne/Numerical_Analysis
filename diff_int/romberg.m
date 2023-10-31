% [s,T] = romberg(funct,a,b,tol,N)
%
%  We use Romberg Integration combined with the trapezoid rule
%  to approximate the integral of a function.
%
%  The following must be given:
%    The function given in  "funct",
%    The endpoints  a  and  b  of the interval  [a,b],
%    The tolerance tol for the acceptable error; namely,
%      the iteration will stop when |T(k-1,k-1) - T(k,k)| < tol .
%    The maximum number  N  of iterations before giving up on
%      approximating the integral.
%
%  The program gives an approximation  s  of the integral.
%
function [s,T] = romberg(funct,a,b,tol,N)
  % default arguments
  arguments
    funct;
    a;
    b;
    tol double = 0.0001;
    N = 20;
  end

  s = NaN;
  flag = 0;
  T(1,1) = trapezoid(funct,a,b,1);
  
  for n=2:1:N
      T(n,1) = trapezoid(funct,a,b,2^(n-1));
      for j=2:1:n
          T(n,2*j-1) = ((4^(j-1))*T(n,2*j-3) - T(n-1,2*j-3))/(4^(j-1)-1);
          if ( n -j > 0 )
              T(n-1,2*j-2) = (T(n-1,2*j-3) - T(n-2,2*j-3))/(T(n,2*j-3)- T(n-1,2*j-3));
          end
      end
      if ( abs( T(n,2*n-1) - T(n-1,2*n-3) ) < tol )
          flag = 1;
          break;
      end
  end

  if ( flag == 1 )
      s = T(n,2*n-1);
  else
      disp('Cannot satisfy the required toleracne for the error.')
  end
end
