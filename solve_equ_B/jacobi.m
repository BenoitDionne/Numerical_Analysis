%   xx = jacobi(A,b,x,tol,limit)
%
%   We use Jacobi Iterative Method, to find a solution of  A x = b.
%
%   The following must be given:
%     The (n x n)-matrix A.
%     The n-dimensional column vector  b.
%     The n-dimensional column vector  x  ; the first guess.
%     The tolerance  tol .
%     The maximal number  limit  of iterations before to quit.
%
%   The program stops as soon as one of the following conditions
%   is satisfied:
%       The total number of iterations is greater than  limit  or
%       the difference between two successive iterations is
%       less than the maximal tolerated error  tol .
%
function xx = jacobi(A,b,x,tol,limit)
  xx = NaN;
  dim = size(A,1);

  for k = 1:dim
    if ( A(k,k) == 0 )
      disp 'The Jacobi Iterative Method fails because some of the elements on the diagonal are zero.'
      return;
    end
  end

  for k = 1:limit
    xx(1,1) = (b(1,1) - A(1,2:dim)*x(2:dim,1))/A(1,1);
    if dim > 2
      for m = 2:dim-1
        xx(m,1) = (b(m,1) - A(m,1:m-1)*x(1:m-1,1) - A(m,m+1:dim)*x(m+1:dim,1))/A(m,m);
      end
    end
    xx(dim,1) = (b(dim,1) - A(dim,1:dim-1)*x(1:dim-1,1))/A(dim,dim);

    if ( norm(xx - x) < tol)
      disp(sprintf('Number of iterations = %d',k))
      return;
    end

    x=xx;
  end

  disp 'The Jacobi Iterative Method failed to give an approximation to a solution of A x = b within the required accuracy.'
  xx = NaN;
end
