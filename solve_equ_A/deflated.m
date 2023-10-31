%   b = deflated(a,x0)
%
%   We use Horner's method to compute the deflated polynomial
%   q in f(x) = (x-x0) * q(x) .  x0 must be a root of the
%   polynomial f.  f  is a polynomial of degree  d .
%
%   The following must be given:
%     The coefficients  a(i)  of the polynomial
%       f(x) = a(1) + a(2) * x + a(3) * x^2 + ... + a(d+1) * x^d,
%     The root x0,
%
%   The program will compute the coefficients  b(1) , b(2) , ... ,
%      b(d)  of the deflated polynomial
%      q(x) = b(1) + b(2) * x + ... + b(d) * x^(d-1)
%   The program will return b.
%
function b = deflated(a,x0)
  d = length(a)-1;
  b(d) = a(d+1);

  for i=d-1:-1:1
    b(i) = a(i+1) + x0 * b(i+1);
  end
end
