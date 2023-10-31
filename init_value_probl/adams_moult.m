%  [t,w] = adams_moult('function',t0,y0,tf,N)
%
%  This program uses Adams-Moulton predictor-corrector method
%  to solve a ODE
%     y' = f(t,y)
%  with initial condition y(t_0) = y_0.
%
%  The following must be given:
%    The initial conditions  t0  and  y0,
%    The final value  tf  of  t
%    The number  N  of subintervals of  [t0,tf] .
%    The function  f(t,y)  must be given in  funct.
%
%  The program gives:
%    t(i) and the apxroximation  w(:,i)  of  y( t(i) ) 
%    where  t(i) = t0 + h*i  with  h=(tf-t0)/N  and i = 0,1,...,N .
%
function [t,w] = adams_moult(funct,t0,y0,tf,N)
    h = (tf - t0 )/ N;
    t(1) = t0;
    w(1) = y0;

    [t,w] = rgkt4(funct,h,3,t(1),w(1));

    for i=1:4
        ww(i) = w(i);
        f(i) = feval(funct,t(i),ww(i));
    end
    
    for i=4:N
        ttt =  t0 + i*h;
        predict= ww(4) + h*(55.0*f(4) - 59.0*f(3) + 37.0*f(2) - 9.0*f(1))/24.0;
        f(5) = feval(funct,ttt,predict);
        correct = ww(4) + h*(9.0*f(5) + 19.0*f(4) - 5.0*f(3) + f(2))/24.0;
        ww(5) = correct;
        f(5) = feval(funct,ttt,correct);
        ww(1:4) = ww(2:5);
        f(1:4) = f(2:5);
        t(i+1) = ttt;
        w(i+1) = ww(4);
    end
end
