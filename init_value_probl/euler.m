%  [t,w] = euler(funct,t0,y0,tf,npoints)
%
%  We use  Euler's algorithm to solve an  ODE  of the form
%    y = f(t,y)
%  with initial conditons  y(t0) = y0 .
%
%  The following must be given:
%    The initial conditions   y0 = w0  at  t0.
%    The value  tf  of  t  where the integration ends
%    The number  N  of subintervals of the interval  [t0,tf]
%    The function  f(t,y)  must be given in  funct.
%
%  The program will produce an approximation  w_i  of  y(t0+i*h)
%  for  i = 0 to N --- h=(tf-t0)/N .
%
function [t,w] = euler(funct,t0,w0,tf,N)
    h = (tf-t0)/N;
    t(1)=t0;
    w(1)=w0;
    for i=1:N
        w(i+1) = w(i) + h * funct(t(i),w(i));
        t(i+1) = t0 + i*h;
    end
end
