%  function adams_bashforth_4s(funct, sol, t0, tf, N, x0)
%    funct is the right hand side  in y' = f(t,y)
%    sol is the solution for the initial condition y(t_0) = y_0
%    t0 is t_0
%    tf is t_f
%    N is the number of intervals
%    y0 is y_0
%
% Example:
% ff = @(t,y) 1 - 2*y;
% t0=0; t1=4; N=128; y0 = 1;
% [tt,ww] = adams_bashforth_4s(ff, t0, t1, N, y0)
%
function [tt,ww] = adams_bashforth_4s(funct, t0, t1, N, y0)
    h = (t1-t0)/N;
    [tt,ww] = rgkt4(funct,h,3,t0,y0);
    fim3 = funct(tt(1),ww(1));
    fim2 = funct(tt(2),ww(2));
    fim1 = funct(tt(3),ww(3));
    fi = funct(tt(4),ww(4));
    for j=4:1:N;
        index = j+1;
        tt(index) = tt(1)+j*h;
        ww(index) = ww(j) + h*(55*fi - 59*fim1 + 37*fim2 - 9*fim3)./24;
        fim3 = fim2;
        fim2 = fim1;
        fim1 = fi;
        fi = funct(tt(index),ww(index));
    end
end
