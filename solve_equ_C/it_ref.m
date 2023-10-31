%  x = it_ref(A,b,tol,limit)
%
%  We use gaussian elimination with scaled column pivoting and
%  iterative refinement to solve a system of linear equations of
%  the form
%    A(1,1)*x(1,1)   + ... + A(1,dim)*x(dim,1)   = b(1)
%    . . .
%    A(dim,1)*x(1,1) + ... + A(dim,dim)*x(dim,1) = b(dim)
%
%  The following must be given:
%     The  (dim x dim) matrix A,
%     The column vector  b ,
%     The tolerance  tol for the accuracy,
%     The maximal number  limit  of iterations in the
%       iterative refinement.
%
%  The program gives an approximation  x(i,1) , for i =1,
%  ... , dim,  of the solution of this system if the solution
%  is unique.  The program uses iterative refinement to improve the
%  approximation of the solution given by gaussian elimination.
%
%  Example:
%  A = [[0.04,0.01,-0.01]; [0.2, 0.5, -0.2];[1,2,4]]
%  b = [0.0601 ; 0.302 ; 11.03]
%  it_ref(A,b,10^(-7),100)
%
function x = it_ref(A,b,tol,limit)
    dim = size(A,1);
    x = NaN;

    % We save  A  and  b  for the iterative rafinement.
    bb = b;
    AA = A;

    for i=1:dim
        nrow(i) = i;
    end
    
    % We use gaussian elimination to write the system in echelon form.
    for j=1:(dim-1)
        order = j;
        rowmax = norm(A(nrow(j),j:dim), inf);
        if (rowmax == 0)
            disp 'The matrix  A  is not invertible.';
            return;
        end
        max = abs( A(nrow(j),j) )/rowmax;
        for i=(j+1):dim
            rowmax = norm(A(nrow(i),j:dim), inf);
            if (rowmax == 0)
                disp 'The matrix  A  is not invertible.';
                return;
            end
            test = abs( A(nrow(i),j) )/rowmax;
            if (test > max)
                max = test;
                order = i;
            end
        end

        if (max == 0)
            disp 'There is no unique solution.';
            return;
        end

        if (j ~= order)
            ncopy = nrow(j);
            nrow(j) = nrow(order);
            nrow(order) = ncopy;
        end

        for i=(j+1):dim
            A(nrow(i),j) = A(nrow(i),j)/A(nrow(j),j);
            A(nrow(i),(j+1):dim)=A(nrow(i),(j+1):dim) - A(nrow(i),j)*A(nrow(j),(j+1):dim);
            b(nrow(i))=b(nrow(i)) - A(nrow(i),j)*b(nrow(j));
        end
    end

    if (A(nrow(dim),dim) == 0)
        disp 'The matrix  A  is not invertible.';
        return;
    end

    % We now use backward substitution to get an approximation of the
    % solution of the system.

    x(dim,1) = b(nrow(dim))/A(nrow(dim),dim);
    for i=(dim-1):-1:1
        x(i,1) = b(nrow(i));
        x(i,1) = x(i,1) - A(nrow(i),(i+1):dim)*x((i+1):dim,1);
        x(i,1) = x(i,1)/A(nrow(i),i);
    end

    %  We refine the solution of  A x = b  using iterative refinement.

    for (n=1:limit)
        r = bb - AA*x;

        % We use gaussian elimination to solve  A xx = r .
        for (j=1:dim)
            r(nrow(j+1:dim)) = r(nrow(j+1:dim)) - A(nrow(j+1:dim),j)*r(nrow(j));
        end

        % We now use backward substitution to get an approximation of the
        % solution of the system  A xx = r.
        xx(dim,1) = r(nrow(dim))/A(nrow(dim),dim);
        for i=(dim-1):-1:1
            xx(i,1) = r(nrow(i));
            xx(i,1) = xx(i,1) - A(nrow(i),(i+1):dim)*xx((i+1):dim);
            xx(i,1) = xx(i,1)/A(nrow(i),i);
        end
        x = x + xx;

        % We determine if we should stop the iterative refinement.
        if ( norm(xx,inf) < tol )
            return;
        end
    end

    x = NaN;
    disp(sprintf('The iterative refinement failed to give an approximation ',...
                 'to a solution of  A x = b  after %d iterations',limit))
end
