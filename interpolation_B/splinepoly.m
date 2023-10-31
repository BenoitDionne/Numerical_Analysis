%  [y, coeffs] = splinepoly(z,f,x,X)
%
%  This programme computes the cubic spline poynomial
%    defined by the solution of the system A z = b where A
%    is the tridiagonal matrix given by splinematrix(), and
%    f(i) = f(x(i)) for i = 1 to npoints.
%
%  The following must be given:
%    The points x(i) for i=1 to npoints.
%    The values f(i) = f(x(i)) for i = 1 to npoints.
%    The solution z of the system Az = b associated to the cubic spline used.
%    The values of x where the cubic spline interpolation polynomial must
%      be evaluated (X in the code below).
%
%  This function returns:
%    The value of the cubic spline polynomial at all the given values of x.
%    The coefficients for each polynomials
%      p_i(x) = ((coeffs(i,1) (x-x(i)) + coeffs(i,2))(x-x(i))
%             + coeffs(i,3))(x-x(i)) + coeffs(i,4)
%      for i=1,2,...,npoins-1 .
%
function [y, coeffs] = splinepoly(z,f,x,X)
  npoints = length(x);
  N = length(X);
  y = repmat(NaN,1,N);

  for m=1:1:npoints-1
    coeffs(m,4) = f(m);
    dx = x(m+1)-x(m);
    df = f(m+1)-f(m);
    coeffs(m,3) = -(2*z(m)+z(m+1))*dx/6 + df/dx;
    coeffs(m,2) = z(m)/2;
    coeffs(m,1) = (z(m+1)-z(m))/(6*dx);
  end

  for n=1:1:N
    J = 0;
    if ( X(n) >= x(1) && X(n) <= x(npoints) )
      for m = 2:1:npoints
        if ( X(n) <= x(m) )
          J = m-1;
          break;
        end
      end
      dx = X(n) - x(J);
      y(n) = ((coeffs(J,1)*dx + coeffs(J,2))*dx + coeffs(J,3))*dx + coeffs(J,4);
    end
  end
end
