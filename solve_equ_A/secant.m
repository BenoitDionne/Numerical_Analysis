%  z = secant(funct,x,y,tol,limit)
%  We use the Secant Method to approximate a root of the equation f(x) = 0.
%
%  The following must be given:
%    The function  f  is given in  funct.
%    The two first guesses  x  and  y.
%    The tolerance  tol
%    The maximal number  limit  of iterations.
%
%   The program stops as soon as one of the following conditions
%   is satisfied:
%     the total number of iterations is greater than limit.
%     the difference between two iterations is less than the tolerance  tol.
%     |f(x) - f(y)| is smaller than the tolerance  tol.
%
% Example:
%  f = @(x) exp(x) - tan(x);
%  x = 1.3;
%  y = 1.35;
%  tol = 10^(-8);
%  limit = 100;
%  secant(f,x,y,tol,limit)
%

function z = secant(funct,x,y,tol,limit)
  z = NaN;
  v = funct(x);
  p = y - x;
  disp(sprintf('n= %d   x= %.12f   y= %.12f   |y-x|= %.12f',0,x,y,p)); 

  for n = 1:limit
    if (abs(p) < tol)
      disp 'We have that |x - y| < tol.'
      z = y;
      return;
    end

    w = funct(y);
    q = w - v;
    if ( abs(q) < tol )
      disp 'We have that  f(x)  is almost equal to  f(y)'
      disp '(i.e.  |f(x) - f(y)| < tol), thus we are close'
      disp 'to a division by zero.'
      return;
    end

    p = -w*p/q;
    z = y + p;

    disp(sprintf('n= %d   x= %.12f   y= %.12f   |y-x|= %.12f',n,y,z,p)); 
    
    v=w;
    y=z;
  end

  disp(sprintf('After %d iterations, the secant method failed to give an approximation of the root of  f(x).',limit));
end
