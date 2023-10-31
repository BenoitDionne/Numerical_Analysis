%  y = interpolation(d,x)
%
%  The following must be given:
%    The matrix d = [ x1, g1 ; x2, g2 ; ...]
%      where  x1 <= x2 <= x3 ...  and
%      g1 = f(x1),
%      g2 = f(x2) if x1 < x2  or  f'(x2) if x1 = x2 
%      g3 = f(x3) if x2 < x3, f'(x3) if x1 < x2 = x3  or  f''(x3) if
%        x1 = x2 = x3
%      In general gk = f^(j)(xk) if xk is repeated j times before.
%    The value x where to evaluate the interpolating polynomial
%
%  The function returns:
%    The value of the interplating polynomial at  x
%
function y = interpolation(d,x)
  n = size(d,1);
  [z,table] = divideddiff(d);
  y = polynomial(d(1:n-1,1)',z,x);         
end
