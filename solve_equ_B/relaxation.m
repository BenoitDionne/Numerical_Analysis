%   xx = relaxation(A,b,x,w,tol,limit)
%
%   We use the relaxation method, to find a solution of  A x = b.
%
%   The following must be given:
%     The (n x n)-matrix A.
%     The n-dimensional column vector  b.
%     The n-dimensional column vector  x  ; the first guess.
%     The number w.
%     The tolerance  tol .
%     The maximal number  limit  of iterations before to quit.
%
%   The program stops as soon as one of the following conditions
%     is satisfied:
%     The total number of iterations is greater than limit  or
%     the difference between two iterations is less than
%     the maximal tolerated error  tol .
%
function xx = relaxation(A,b,x,w,tol,limit)
  xx = NaN;
  dim = size(A,1);

  for k = 1:dim
    if ( A(k,k) == 0 )
      disp 'The Relaxation Iterative Method fails because some of the'
      disp 'elements on the diagonal are zero.'
      return;
    end
  end

  for k = 1:limit
    xx(1,1) = x(1,1) + w*(b(1,1) - A(1,:)*x(:,1))/A(1,1);
    if dim > 2
      for m = 2:dim
        xx(m,1) = x(m,1) + w*(b(m,1) - A(m,1:m-1)*xx(1:m-1) - A(m,m:dim)*x(m:dim))/A(m,m);
      end
    end

    if ( norm(xx - x) < tol)
        disp(sprintf('Number of iterations = %d',k))
        return;
    end
    
    x=xx;
  end

  disp 'The Jacobi Iterative Method failed to give an approximation to a'
  disp 'solution of  A x = b  within the rquired accuracy.'
  xx = NaN;
end
