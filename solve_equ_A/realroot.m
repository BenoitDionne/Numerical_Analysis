%   xf = realroot(a,x,N,tol)
%
%   We use Newton-Raphson to approximate a real root of a
%   polynomial p(x) = a_n x^n + ... a_1 x + a_0 using Horner's algorithm
%
%   Input:
%     The coefficients a_i (the vector a  in the code below),
%       (even the zero coefficients must be given.)
%     The initial approximation x_0 (x in the code below) of a root
%       c of p
%     The maximal tolerance  tol
%     The maximal number N of iterations
%   Output:
%     An approximation (xf in the code below) of the
%     real root c  or  an error message if the real root cannot be
%     approximate with the desired tolerance in less than N iterations.
%
function xf = realroot(a,x,N,tol)
  xf = NaN;
  m = length(a);

  for k=1:N
    y = a(m);
    z = a(m);
    for i=m-1:-1:2
      y = a(i) + x*y;
      z = y + x*z;
    end
    y = a(1) + x*y;

    if ( abs(z) < tol )
      disp 'The derivative is almost null.  Cannot proceed.'
      return;
    end

    % y = f(x)  and  z = f'(x) .
    ratio = y/z;
    x = x - ratio;
    if (abs(ratio) < tol)
      xf = x;
      return;
    end
  end
  
  disp 'The program fails to give an approximation to a root of'
  disp 'the polynomial in less than the maximum number of iterations.'
  xf = NaN;
end
