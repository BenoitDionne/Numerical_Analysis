% E = QRshifting(a, b, N, d)
%
%   We use the QR algorithm with shisting to approximate the
%   eigenvalues of a tridiagonal symmetric matrix .
%
%   The following must be given:
%     The vector of elements on the diagonal  a,
%     The vector of element on the subdiagonal b,
%     The maximal number of iteration for the QR algorithm
%     The function  f  is given in  funct.
%
%   The program return:
%     The approximation of the eigenvalues.

function E = QRshifting(a, b, N, d)
  n = length(a);
  E = repmat(NaN,n,1);
  nE = 0;
  s = 0;          % sum of the shifts

  for k = 1:N
    fprintf('%d ... ',k);

    % If n=1, we are done
    if ( n == 1 )
      nE = nE + 1;
      E(nE) = a(1) + s;
      return;
    end

    % If the matrix can be splitted into two symmetric tridiagonal
    % matrices, we do so.
    for j = 1:n-1
      if ( abs(b(j)) < d )
        disp 'Splitting the matrix';
        E(nE+1:nE+j,1) = QRshifting(a(1:j),b(1:j-1), N, d) + s;
        E(nE+j+1:n,1) = QRshifting(a(j+1:n),b(j+1:n-1), N, d) + s;
        return;
      end
    end

    % We compute the eigenvalues of the matrix
    %    [ a_{n-1} b_{n-1} ]
    %    [ b_{n-1} a_n     ]
    % We use the appropriate form of the formula to find the roots
    % of a quadratic equation to avoid subtraction of almost equal
    % numbers.
    B = -(a(n-1) + a(n));
    C = a(n)*a(n-1) - b(n-1)*b(n-1);
    D = sqrt(B^2-4*C);
    if ( B > 0 )
      r1 = -2*C/(B+D);
      r2 = -(B+D)/2;
    else
      r1 = (D-B)/2;
      r2 = 2*C/(D-B);
    end

    % If we have only a 2 x 2 matrix, we have found approximations
    % for the last two eigenvalues of A.
    if ( n == 2)
      nE = nE + 1;
      E(nE,1) = r1 + s;
      nE = nE + 1;
      E(nE,1) = r2 + s;
      return;
    end

    % Chose the appropriate shift
    if ( abs(r1-a(n)) < abs(r2 -a(n)) )
      stmp = r1;
    else
      stmp = r2;
    end
    s = s + stmp;
    a = a - stmp;

    % Get the QR decomposition
    x = a(1);
    y = b(1);
    for j = 1:n-1
      alpha(j) = sqrt( x^2 +b(j)^2 );
      ccc(j) = x/alpha(j);
      sss(j) = b(j)/alpha(j);
      beta(j) = y*ccc(j) + a(j+1)*sss(j);
      x = - y*sss(j) + a(j+1)*ccc(j);
      if ( j ~= n-1 )
        gamma(j) = b(j+1)*sss(j);
        y = b(j+1)*ccc(j);
      end
    end
    alpha(n) = x;

    % Compute RQ knowing that the result is a symmetric tridiagonal matrix
    a(1) = alpha(1)*ccc(1) + beta(1)*sss(1);
    b(1) = alpha(2)*sss(1);
    for j = 2:n-1;
      a(j) = alpha(j)*ccc(j-1)*ccc(j) + beta(j)*sss(j);
      b(j) = alpha(j+1)*sss(j);
    end
    a(n) = alpha(n)*ccc(n-1);
  end
end
