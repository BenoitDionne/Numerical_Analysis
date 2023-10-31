%  z = tridmatrix(L,D,U,b)
%
%  We solve a tridiagonal system of the form  A x = b.
%
%  The function must be given:
%    The lower diagonal L, the diagonal D and the upper diagonal U of
%      the tridiagonal matrix  A.  None of the components of
%      the diagonal $D$ can be null.
%    The right hand side b
%
%  The function returns:
%     The solution if the system can be solved.
%
function [z,diagnostic] = tridmatrix(L,D,U,b)
  m = length(D);
  z = repmat(NaN,1,m);

  for n=2:m
    if (D(n-1) == 0)
      return;
    end
    q = L(n-1)/D(n-1);
    D(n) = D(n)-q*U(n-1);
    b(n) = b(n)-q*b(n-1);
  end

  if (D(m) == 0)
    return;
  end
  
  % Backward substitution
  z(m) = b(m)/D(m);
  for n=(m-1):-1:1
    z(n)=(b(n)-U(n)*z(n+1))/D(n);
  end
end
