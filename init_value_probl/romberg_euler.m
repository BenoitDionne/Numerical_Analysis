%  [gt,gw] = ronmberg_euler(funct,t0,w0,tf,hmin,hmax,tol,N)
%
%  This program uses a Romberg Integration between two
%  consecutive nodes to approximate solution of the ODE
%    y' = f(t,y)
%  with initial condition y(t0) = y0.  Euler Method is used
%  generate the first column of the Romberg table.  Iterations
%  for the Romberg table end when |T(n,n) - T(n-1,n-1)| < tol.
%  The value T(n,n) is used as the approximation w(i+1) of y(t(i+1)).
%       
%  Since the local discretization error of the Euler Method
%  is not a sum of terms which are factor of h^(2k), we
%  don't have division by a factor of h^2 of the error but only
%  a factor of h.
%
%  The following must be given:
%    The maximal step size  hmax,
%    The minimal step size  hmin,
%    The maximal tolerated error  tol,
%    The initial conditions  t0  and  y0 = w0,
%    The final value  tf  of  t  where to stop the integration.
%    The function  f(t,y)  must be given in  funct.
%    The maximal number N of iterations for the Romberg table 
%
%  There are two flags:
%    flag = 1  if Romberg Integration Method is successful and
%    flag = 0  otherwise.
%    last = 1  when we reach  tf  or  h < hmin, and last = 0  otherwise.
%
%  The program gives:
%    gt(i) and the apxroximation gw(i) of y(t(i)) .
%
function [gt,gw] = romberg_euler(funct,t0,w0,tf,hmin,hmax,tol,N)
    % default arguments
    arguments
        funct;
        t0;
        w0 double;
        tf;
        hmin double = 0.001;
        hmax double = 0.1;
        tol double = 0.0001;
        N = 20;
    end

    k = 1;
    gt(k) = t0;
    gw(k) = w0;
    last = 0;
    h=hmax;

    while (last == 0)
        flag = 0;
        t1 = gt(k);
        t2 = gt(k)+h;
        [tt,ww] = euler(funct,t1,gw(k),t2,1);
        T(1,1) = ww(length(ww));

        for n=2:1:N
            [tt,ww] = euler(funct,t1,gw(k),t2,2^(n-1));
            T(n,1) = ww(length(ww));
            for j=2:1:n
                T(n,j) = ((4^(j-1))*T(n,j-1) - T(n-1,j-1))/(4^(j-1)-1);
            end
            if ( abs( T(n,n) - T(n-1,n-1) ) < tol )
                flag = 1;
                break;
            end
        end

        if ( flag == 1 )
            k = k+1;
            gt(k) = t2;
            gw(k) = T(n,n);
            if ( abs(tf-t2) < hmin/2)
                last = 1;
            elseif (t2 + h > tf)
                h = tf - t0;
            elseif ( (n<N) & (h < hmax/2) & (t2 + 2*h <= tf) )
               h = 2*h;
            end
        else
            h = h/2.0;
            if (h < hmin)
                disp 'The step size has to be smaller than hmin.'; 
                gt=NaN;
                gw=NaN;
                return;
            end
        end
    end
end
