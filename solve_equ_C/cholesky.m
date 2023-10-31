%  x = cholesky(A,b)
%
%  We use Cholesky Factorzation to solve a system of linear
%  equations of the form
%
%    A(1,1)*x(1)     + ... + A(1,dim)*x(dim)   = b(1,:)
%    . . .
%    A(dim,1)*x(x(1) + ... + A(dim,dim)*x(dim) = b(dim,:)
%
%  The following must be given:
%    The square matrix A
%    The matrix ( b(:,i) ) for i=1, 2, ..., M ; the M linear
%      systems A x = b(:,i) are solved simultaneously.
%
%  The program gives an approximation  x(:,i)  of the solution of
%  the linear system associated to b(:,i) for i=1, 2, ..., M.
%
function [x,M] = cholesky(A,b)
  dim = size(A,1);

  % In theory, we do not have to use pivoting.  Moreover, we only need
  % to compute  M  because  N  is the trampose of  M.
  M = zeros(dim,dim);
  
  M(1,1) = sqrt(A(1,1));
  M(2:dim,1) = A(2:dim,1)/M(1,1);
  for k = 2:dim-1
    M(k,k) = sqrt(A(k,k) - sum(M(k,1:k-1).^2));
    for j = k+1:dim
      M(j,k) = (A(j,k) - sum(M(j,1:k-1).*M(k,1:k-1)))/M(k,k);
    end
  end
  M(dim,dim) = sqrt(A(dim,dim) - sum(M(dim,1:dim-1).^2));

  % Only at this point do we need the value of  b .
  % We use forward substitution to sole  My = b.

  y(1,1) = b(1,1)/M(1,1);
  for i = 2:dim
    y(i,1) = (b(i,1) - M(i,1:i-1)*y(1:i-1,1) )/M(i,i);
  end

  % We now use backward substitution to get an approximation of the
  % solution of the system.

  x(dim,1) = y(dim,1)/M(dim,dim);
  for i = dim-1:-1:1
    x(i,1) = (y(i,1)-M(i+1:dim,i)'*x(i+1:dim,1))/M(i,i);
  end
end
