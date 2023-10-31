%  [t,w] = trapez1(funct,Dfunct,t0,y0,tf,N,max,tol)
%
%  We use the implicit trapezoid method to approximate the
%  solution of the stiff differential equation  y' = f(t,y)
%  with the initial conditions  y(t_0) = y_0 .
%
%  The following must be given:
%   - The initial time t_0 (t0 in the code below)
%   - The final time t_f (tf in the code below)
%   - The number N of subintervals of the interval [t_0,t_f]
%   - The initial conditions y_0 at t_0 (y0 in the code below)
%   - The maximal tolerance tol that we are willing to accept,
%   - The maximal number of iterations max allowed for the
%     Newton's method.
%   - The function f(t,y) (funct in the code below)
%   - The Derivative of f(t,y) with respect to y (Dfunct
%     in the code below)
%
%  The program gives:
%    t(i) = t_0 + i*h and the approximation  w_i  of y_i.
%
%  NOT FULLY TESTED
function [t,w] = trapez1(funct,Dfunct,t0,y0,tf,N,max,tol)
    h = (tf-t0)/N;
    half = h/2;
    t(1) = t0;
    w(:,1) = y0;
    n = length(y0);
    
    for i=1:N
        % We start the iteration with the approximation of y(t_0 + h)
        % given by the Euler method.
        k = feval(funct,t(i),w(i));
        w0 = w(i) + h.*k;
        t(i+1) = t0 + i*h;

        for (j=1:max)
            Fw0 = w0 - w(i) -half*(funct(t(i+1),w0) + k);
            DFw0 = eye(n) - half*Dfunct(t(i+1),w0);
            if (rank(DFw0) ~= n)
                disp 'Newton-Raphson iterative method does not work.'
                t = NaN;
                w = NaN;
                return;
            end
            ww = linsolve(DFw0,Fwo);
            w1 = w0 - ww;
            if (norm(w1-w0) < tol)
                w(:,i+1) = w1;
                break;
            else
                w0 = w1;
                if (j == max)
                    disp 'The maximum number of iterations has been reached';
                    disp 'before being anble to satisfy the tolerance required.';
                    t = NaN;
                    w = NaN;
                    return;
                end
            end
        end
    end
end

