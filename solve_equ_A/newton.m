%  xx = newton(funct,functprime,x,tol,limit)
%  We use Newton's Method to approximate a root of the
%  equation f(x) = 0.
%
%  The following must be given:
%    The function  f  in  funct.
%    The function  f'  in  functprime.
%    The initial guess  x  for a root of  f.
%    The tolerance  tol.
%    The maximum number  limit  of iterations.
%
%  The program stops as soon as one of the following conditions
%  is satisfied:
%    the total number of iterations is greater than limit
%    the difference between two iterations is less than the tolerance  tol .
%    |f'(z)| < tol
%
% Example:
%  f = @(x) atan(x) - 2*x./(1+x.^2);
%  fp = @(x) (-1 + 3*x.^2)/((1+x.^2).^2);
%  x0 = 1;
%  tol = 10^(-12);
%  limit = 100;
%  xx = newton(f,fp,x0,tol,limit)
%
% Example:
%  f = @(x) tan(x);
%  fp = @(x) (sec(x)).^2;
%  x = 5;
%  tol = 10^(-8);
%  limit = 100;
%  newton(f,fp,x,tol, limit)
%
function xx = newton(funct,functprime,x,tol,limit)
  for n=1:limit
    z = functprime(x);
    
    % Instead of testing if  f'(x) = 0, we test if
    % abs(fx) < 2*realmin , where realmin is the smallest number
    % that the computer may handle
    if ( abs(z) < 2*realmin )
      disp 'The derivative is almost zero at some point';
      xx = NaN;
      return
    end

    y = x - funct(x)/z;

    w = abs(x-y);
    disp(sprintf('n= %d   x= %.12f   y= %.12f   |y-x|= %.12f',n,x,y,w)); 

    % if  |x-y| < tol, we may assume that up to the requested
    % accuracy, we have reached a root of f.
    if (w < tol)
      xx = y;
      disp(sprintf('%d iterations were needed.',n));
      return
    end
    
    x=y;
  end
  
  disp(sprintf('Newton-Raphson method failed to give an approximation of a root of  f(x)  after %d iterations.',limit));
  xx = NaN;
end

