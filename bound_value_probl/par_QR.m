% [Q,R] = par_QR(A)
%
%   Find the QR decomposition of a matrix A; namely, A = Q R, where the columns of
%   Q are orthonormal and  R is triangular superior.   The columns of
%   A must be linearly independent. 
%
%   The following must be given:
%     The n times q matrix A.
%
%   The program return:
%     The Q R decomposition
%
function [Q,R] = par_QR(A)
  n = size(A,1);
  q = size(A,2);
  Q = repmat(NaN,n,q);
  R = zeros(q);
  
  if ( rank(A) < q )
      disp 'The columns of A must be linearly independent.';
      return;
  end

  R(1,1) = norm(A(:,1));
  Q(:,1) = (1/R(1,1))*A(:,1);
  for i = 2:1:q
      v = A(:,i);
      for j = 1:1:i-1
          R(j,i) = Q(:,j)'*A(:,i);
          v = v - R(j,i)*Q(:,j); 
      end
      R(i,i) = norm(v);
      Q(:,i) = (1/R(i,i))*v;
  end
end
