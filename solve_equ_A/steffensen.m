%  xx = steffensen(funct,x0,tol,limit)
%
%  We use Steffensen's Method to approximate a fixed point
%  of g.
%
%  The following must be given:
%    The initial approximation  x,
%    The maximal tolerated error  tol,
%    The maximal number  limit  of iterations.
%    The functions  g  must be given in  funct.
%
%  The program gives an approximation   xx  of a fixed point
%  of g or an error message.
%
%  The program stops as soon as one of the following conditions
%  is satisfied:
%  The total number of iterations is greater than limit or
%  the difference between two iterations is less than the maximal
%  tolerated error tol .
%
% Example:
%  f = @(x) 12 -20./x;
%  x = 9.5;
%  tol = 10^(-7);
%  limit = 100;
%  p = steffensen(f,x,tol,limit)
%
function xx = steffensen(funct,x,tol,limit)
    list = x;      % a list of iterations

    for n=1:limit
        y = feval(funct,x);
        z = feval(funct,y);
        u = y - x;
        w = z - y - u;
        
        if ( abs(w) < 2*realmin )
            disp 'The program has to stop because the forward';
            disp 'difference of power two for  x  is almost zero.';
            xx = x;
            list
            return;
        end

        u = x - u.^2./w;
        list = [ list ; u ];

        if (abs(x-u) < tol)
            xx = u;
            list
            return;
        end

        x=u;
    end
  
    disp(sprintf('After %d iterations, Steffensen''s method failed to give an approximation to a fixed point of g within the required tolerated error.',limit))
    xx = NaN;
end
