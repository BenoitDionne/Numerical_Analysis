%  [t,w] = midpoint(funct,t0,w0,tf,npoints)
%
%  We use midpoint method to solve an  ODE  of the form
%    y = f(t,y)
%  with initial conditons  y(t_0) = y_0
%
%  The following must be given:
%    The initial conditions   y0 = w0  at  t0.
%    The value  tf  of  t
%    The number  N  of subintervals of the interval  [t0,tf]
%    The function  f(t,y)  must be given in  funct.
%    
%  The program will produce an approximation  w(:,i)  of  y(t0+i*h)
%       where  i = 0 to N  and  h=(tf-t0)/N .
%
function [t,w] = midpoint(funct,t0,y0,tf,N)
    h = (tf-t0)/N;
    t(1) = t0;
    w(:,1) = y0;

    for i=1:N
        tempo = w(:,i) + (h/2)*funct(t(i),w(:,i));
        w = [w,w(:,i) + h * funct(t(i)+(h/2),tempo)];
        t(i+1) = t0 + i*h;
    end
end

