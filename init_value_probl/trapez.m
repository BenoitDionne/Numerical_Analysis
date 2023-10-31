%  [t,w] = trapez(funct,functprime,t0,y0,tf,n,max,tol)
%
%  We use the implicit trapezoid method to approximate the
%  solution of the stiff differential equation  y' = f(t,y)
%  with the initial conditions  y(t_0) = y_0 .
%
%  The following must be given:
%   - The initial time t_0 (t0 in the code below)
%   - The final time t_f (tf in the code below)
%   - The number n of subintervals of the interval [t_0,t_f]
%   - The initial conditions y_0 at t_0 (y0 in the code below)
%   - The maximal tolerance tol that we are willing to accept,
%   - The maximal number of iterations max allowed for the
%     Newton's method.
%   - The function f(t,y) (funct in the code below)
%   - The Derivative of f(t,y) with respect to y (functprime
%     in the code below)
%
%  The program gives:
%    t(i) = t_0 + i*h and the approximation  w_i  of y_i.
%
% Example:
% ff = @(t,y) 100*y+100*t.^2-2*t-100;
% ffp = @(t,y) 100;
% t0=0; tf=1; y0=1; max=10; tol=0.00001; n=10;
% [t,r] = trapez(ff,ffp,t0,y0,tf,n,max,tol)
% save results.dat [t;r]' -ascii
%
function [t,w] = trapez(funct,functprime,t0,y0,tf,n,max,tol)
    h = (tf-t0)/n;
    half = h/2;
    t(1) = t0;
    w(1) = y0;
    
    for i=1:n
        % We start the iteration with the approximation of y(t_0 + h)
        % given by the Euler method.
        k = feval(funct,t(i),w(i));
        w0 = w(i) + h.*k;
        t(i+1) = t(1) + i*h;

        % Newton-Raphson iterations
        for (j=1:max)
            numer = w0 - w(i) -half*(funct(t(i+1),w0) + k);
            denum = 1.0 - half*functprime(t(i+1),w0);
            if (denum == 0)
                disp 'Newton-Raphson iterative method does not converge (fast enough).'
                t = NaN;
                w = NaN;
                return;
            end
            w1 = w0 - numer/denum;
            if (abs(w1-w0) < tol)
                w(i+1) = w1;
                break;
            else
                w0 = w1;
                if (j == max)
                    disp 'The maximum number of iterations has been reached';
                    disp 'before getting an approximation of w(i+1) within the';
		    disp 'required tolerance.';
                    t = NaN;
                    w = NaN;
                    return;
                end
            end
        end
    end
end
