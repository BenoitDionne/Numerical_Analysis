%  Y = newton_syst(F,DF,X,N,tol)
%
%  Program to compute iterations for the Newton method
%
%  Given:
%    F : a function from R^2 to R^2  (column vectors)
%    DF : the derivative of F
%    X : the initial condition (column vector)
%    N : the maximal number of iterations.  The fefault value is 5.
%    tol : the conditon to end the iteration.  The iterations end
%         when the difference between two successive iterations is less
%         than tol.  The default value is 0.0001 = 10^(-4) .
%
function Y = newton_syst(F,DF,X,N,tol)

  % default arguments
  arguments
    F;
    DF;
    X (:,1) double;
    N int32 = 5;
    tol double = 0.0001;
  end

  for i=1:1:N
    FX = F(X);
    DFX = DF(X);
    Y = X - DFX\FX;
    err = norm(Y-X);
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
