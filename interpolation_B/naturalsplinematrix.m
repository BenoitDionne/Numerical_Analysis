%  [L,D,U,b] = naturalsplinematrix(f,p)
%
%  This programme computes the tridiagonal matrix  A  and
%    the right hand side column vector  b  associated with
%    the natural cubic spline interpolation.
%
%  The following must be given:
%    The points x(i) for i=1 to N.
%    The values f(i) = f(p(i)) for i = 1 to N.
%
%   This function generates:
%     The lower diagonal L, the diagonal D and the upper diagonal U
%       of the tridiagonal matrix associated to the cubic spline.
%     The right hand side b(i) for i=1 to N-2.
%
function [L,D,U,b] = naturalsplinematrix(f,p)
  N = length(p);
  L = repmat(NaN,1,N-3);
  U = repmat(NaN,1,N-3);
  D = repmat(NaN,1,N-2);
  b = repmat(NaN,1,N-2);

  dp = p(2)-p(1);
  if (dp == 0)
    return;
  end
  ratio = (f(2)-f(1))/dp;

  for n=1:N-2
    prevdp = dp;
    dp = p(n+2)-p(n+1);
    if (dp == 0)
      return;
    end
    prevratio = ratio;
    ratio = (f(n+2)-f(n+1))/dp;
    D(n) = 2*(dp+prevdp);
    if ( n < N - 2 )
        U(n) = dp;
        L(n) = dp;
    end
    b(n) = 6*(ratio - prevratio);
  end
end
