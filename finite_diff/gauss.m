%  x = gauss(A,b,option)
%
%  We use gaussian elimination with maximal column pivoting
%  (obtion = 1) or scaled column pivoting (obtion = 2) to solve
%  a system of linear equations of the form
%
%    A(1,1)*x(1)     + ... + A(1,dim)*x(dim)   = b(1,:)
%    . . .
%    A(dim,1)*x(x(1) + ... + A(dim,dim)*x(dim) = b(dim,:)
%
%  The following must be given:
%    The matrix A
%    The matrix ( b(:,i) ) for i=1, 2, ..., M ; the M linear
%      systems A x = b(:,i) are solved simultaneously.
%    The option  option  chosen: option = 1 for partial column
%      pivoting and option = 2 for scaled column pivoting.
%
%  The program gives an approximation  x(:,i)  of the solution of
%  the linear system associated to b(:,i) for i=1, 2, ..., M.
%
function x = gauss(A,b,option)
  dim = size(A,1);
  x = NaN;

  if ( (option ~= 1) & (option ~= 2) )
    disp 'There is no such algorithm.';
    return;
  end

  % To avoid expensive row interchanges, we only interchange the
  % indices of the rows.  We create the vector N = ( 1 2 3 ... dim )
  % to keep track of the permutations of the rows.  N(i) will contain the
  % index of the row in the original matrix  A  which is now located
  % in row  i .
  N = linspace(1,dim,dim);

  % We use gaussian elimination to write the system in echelon form.

  for j=1:(dim-1)
    % If option = 1, then we use the maximal colum pivoting algorithm.
    % If option = 2, then we use the scaled column pivoting algorithm.

    if (option == 1)
      order = j;
      max = abs( A(N(j),j) );
      for i=(j+1):dim
        if (abs( A(N(i),j) ) > max)
          max = abs( A(N(i),j) );
          order = i;
        end
      end
    else
      % We find the index j such that
      % |a^k_{j,k}|/max_{k\leq i \leq n}|a^k_{j,i}| >=
      %   |a^k_{s,k}|/max_{k\leq i \leq n}|a^k_{s,i}|
      % for k <= s <= dim.
      order = j;
      rowmax = norm(A(N(j),j:dim), inf);
      if (rowmax == 0)
        disp 'The matrix  A  is not invertible.';
        return;
      end
      max = abs( A(N(j),j) )/rowmax;
      for i=(j+1):dim
        rowmax = norm(A(N(i),j:dim), inf);
        if (rowmax == 0)
          disp 'The matrix  A  is not invertible.';
          return;
        end
        test = abs( A(N(i),j) )/rowmax;
        if (test > max)
          max = test;
          order = i;
        end
      end
    end

    if (max == 0)
      disp 'There is no unique solution.';
      return;
    end

    % We interchange the k^{th} and j^{th} rows.
    if (j ~= order)
      ncopy = N(j);
      N(j) = N(order);
      N(order) = ncopy;
    end

    for i=(j+1):dim
      A(N(i),j) = A(N(i),j)/A(N(j),j);
      A(N(i),(j+1):dim) = A(N(i),(j+1):dim) ...
                           - A(N(i),j)*A(N(j),(j+1):dim);
      b(N(i),:)=b(N(i),:) - A(N(i),j)*b(N(j),:);
    end
  end

  % We now use backward substitution to get an approximation of the
  % solution of the system.

  if (A(N(dim),dim) == 0)
    disp 'The matrix  A  is not invertible.';
    return;
  end

  x(dim,:) = b(N(dim),:)/A(N(dim),dim);
  for i=(dim-1):-1:1
    x(i,:) = b(N(i),:);
    for j=(i+1):dim
      x(i,:) = x(i,:) - A(N(i),j)*x(j,:);
    end
    x(i,:) = x(i,:)/A(N(i),i);
  end
end
