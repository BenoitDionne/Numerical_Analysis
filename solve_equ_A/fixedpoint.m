%  xx = fixedpoint(funct,x,tol,limit)
%  We use fixed point method, to find a fixed point of a function.
%
%  The following must be given:
%    The initial guess  x  for the fixed point.
%    The tolerance  tol  accepted.
%    The maximum number  limit  of iterations allowed.
%    The function  g  must be given in  funct.
%
%  The program stops as soon as one of the following conditions
%    is satisfied:
%    The total number of iterations is greater than limit or
%    the difference between two consecutive iterations is less than
%    the maximal error  tol .
%
function y = fixedpoint(funct,x,tol,limit)
    list = x;

    for n=1:limit
        y = funct(x);
        list = [ list ; y ];
        if (abs(x-y) < tol)
            list
            disp(sprintf('number of iteration = %d.',n))
            return;
        end
        x = y;
    end

    list
    disp(sprintf('Cannot get an approximation of a fixed point of g after %d iterations',limit))
    y = NaN;
end
