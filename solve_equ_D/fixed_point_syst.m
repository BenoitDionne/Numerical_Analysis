%  Y = fixed_point_syst(G,X,N,tol)
%
%  Program to compute iterations for the Fixed Point method
%
%  Given:
%    G : a function from R^2 to R^2  (column vectors)
%    X : the initial condition (column vector)
%    N : the maximal number of iterations.  The fefault value is 5.
%    tol : the conditon to end the iteration.  The iterations end
%         when the difference between two successive iterations is less
%         than tol.  The default value is 0.0001 = 10^(-4) .
%    D : the norm to use, D=0 for the 2-norm (default), D=1 for the
%        1-norm, and D != 0 or 1 for the infinity norm.
%
%  Example:
% 
%  G = @(x) [(x(1,1).^2+x(2,1).^2+8)/10 ; (x(1,1).*(x(2,1).^2)+x(1,1)+8)/10];
%     Don't leave blank spaces in math formulae used to define functions.
%     Matlab doesn't accept that.
%
%  Y = fixed_point_syst(G, [0;0], 100, 0.00003)
%
function Y = fixed_point_syst(G,X,N,tol,D)

  % default arguments
  arguments
    G;
    X (:,1) double;
    N int32 = 5;
    tol double = 0.0001;
    D int16 = 0;
  end

  for i=1:1:N
    Y = G(X);
    if ( D == 0 )
        err = norm(Y-X);
    elseif ( D == 1 )
        err = norm(Y-X,1);
    else
        err = norm(Y-X,Inf)
    end
    if ( err < tol )
      break
    end
    X = Y;
  end

  disp('X = ')
  X
  disp('Iteration number:')
  i
  disp('Difference from previous iteration:')
  err
end
