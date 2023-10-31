%  v = polynomial(X,a,x)
%
%  Evaluate the polynomial
%    p(x) = a(1) + [a(2) + [a(3) + ...
%        [a(n-1)+a(n)*(x-X(n-1))]*(x-X(n-2))...]*(x-X(2))]*(x-X(1)) .
%    at a given point using Horner Algorithm.
%
%  The following must be given:
%    The coefficients  a(1), a(2), ... a(n)  of the polynomial
%       (the vector a in the code below).
%    The points  X(1),...,X(n-1) of the interpolating polynomial
%       (the vector p in the code below). The point X(n) is not needed.
%    The value  x  where to evaluate the polynomial
%
%  The function returns:
%    The value of the interplating polynomial at  x
%
function v = polynomial(X,a,x)
  n=length(a);
  v = a(n);
  for i=(n-1):-1:1
    v=a(i)+(x-X(i)).*v;
  end
end
