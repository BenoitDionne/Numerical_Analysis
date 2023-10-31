%  function modified_euler(funct, t0, tf, N, x0)
%    funct is the right hand side  in y' = f(t,y)
%    sol is the solution for the initial condition y(t_0) = y_0
%    t0 is t_0
%    tf is t_f
%    N is the number of intervals
%    y0 is y_0
%
% Example:
% format long
% ff = @(t,y) 100*y + 100*t.^2 -2*t -100;
% t0=0; t1=1; N=40; y0 = 1;
% [tt,ww] = modified_euler(ff, t0, t1, N, y0)
%
function [tt,ww] = modified_euler(funct, t0, t1, N, y0)
    h = (t1-t0)/N;
    tt(1) = t0;
    ww(1) = y0;
    for j=1:1:N;
        index = j+1;
        tt(index) = tt(1)+j*h;
        fi = funct(tt(j),ww(j));
        ww(index) = ww(j) + h.*(fi + funct(tt(j)+h, ww(j)+h*fi))/2;
    end
end
