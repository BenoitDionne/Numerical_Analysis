%  [a, table] = divideddiff(d)
%
%  Newton's interpolatory divided-difference formula.
%
%  The following must be given:
%    The matrix d = [ x1, g1 ; x2, g2 ; ...]
%      where  x1 <= x2 <= x3 ...  and
%      g1 = f(x1),
%      g2 = f(x2) if x1 < x2  or  f'(x2) if x1 = x2 
%      g3 = f(x3) if x2 < x3, f'(x3) if x1 < x2 = x3  or  f''(x3) if
%        x1 = x2 = x3
%      In general gk = f^(j)(xk) if xk is repeated j times before.
%
%  The function gives:
%    The table of divided-difference formula
%    The coefficients of the interpolating polynomial
%      p(x)=a(1) + a(2)*(x-x1) + a(3)*(x-x1)*(x-x2) + ... 
%          + a(n)*(x-x1)*...*(x-xm)   with m = n-1.
%
function [a, table] = divideddiff(d)
  n = size(d,1);
  md = 0;
  
  % generate the table of divided differences
  table = repmat(NaN,n,n+1);

  table(:,1) = d(:,1);
  table(1,2) = d(1,2);
  for i=2:n
    if ( table(i,1) == table(i-1,1) )
      table(i,2) = table(i-1,2);
    else
      table(i,2) = d(i,2);
    end
  end
  for k=3:n+1
    for i=1:n-k+2
      m = i+k-2;
      if ( table(m,1) == table(i,1) )
        if ( md == 0 )
          md = m;
        end
        table(i,k) = d(md,2)/factorial(k-2);
      else 
        table(i,k)=(table(i+1,k-1)-table(i,k-1))/(table(m,1)-table(i,1));
        md = 0;
      end
    end
    md = 0;
  end

  a = table(1,2:n+1);
end
