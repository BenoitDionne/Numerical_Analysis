%  [L,D,U,b] = clampedsplinematrix(f,fx,x)
%
%  This programme computes the tridiagonal matrix  A  and
%    the right hand side column vector  b  associated with
%    the clamped cubic spline interpolation.
%
%  The following must be given:
%    The points x(i) for i=1 to N
%    The values f(i) = f(x(i)) for i = 1 to N
%    The values fx(1) = f'(x(1)) and fx(2) = f'(x(N)).
%
%   This function generates:
%     The lower diagonal L, the diagonal D and the upper diagonal U
%       of the tridiagonal matrix associated to the cubic spline.
%     The right hand side b(i) for i=1 to N.
%
function [L,D,U,b] = clampedsplinematrix(f,fx,x)
  N = length(x);
  L = repmat(NaN,1,N-1);
  U = repmat(NaN,1,N-1);
  D = repmat(NaN,1,N);
  b = repmat(NaN,1,N);

  dx = x(2)-x(1);
  if (dx == 0)
    return;
  end
  ratio = (f(2)-f(1))/dx;
  D(1) = 2*dx;
  U(1) = dx;
  b(1) = 6*(ratio - fx(1));

  for n=2:N-1
    prevdx = dx;
    dx = x(n+1)-x(n);
    if (dx == 0)
      return;
    end
    prevratio = ratio;
    ratio = (f(n+1)-f(n))/dx;
    L(n-1) = prevdx;
    D(n) = 2*(dx+prevdx);
    U(n) = dx;
    b(n) = 6*(ratio - prevratio);
  end

  L(N-1) = dx;
  D(N) = 2*dx;
  b(N) = 6*(fx(2) - ratio);
end
