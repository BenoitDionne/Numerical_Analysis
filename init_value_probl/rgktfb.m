%  [gt,gw] = rgktfb(funct,t0,y0,tf,hmin,hmax,T)
%
%  This program uses a Runge-Kutta-Fehlberg method of order four to solve
%  an ODE
%    y' = f(t,y)
%  with initial condition y(t0) = y0.
%  The step size is controlled by a Runge-Kutta method of order five.
%
%  The following must be given:
%    The maximal step size  hmax,
%    The minimal step size  hmin,
%    The maximal tolerated error  T,
%    The initial conditions   y0  at  t0,
%    The final value  tf  of  t  where the integration ends.
%    The function  f(t,y)  must be given in  funct.
%
%  The program gives:
%     t(i) , the apxroximation  w(i)  of  y( t(i) ) ,
%
function [gt,gw] = rgktfb(funct,t0,y0,tf,hmin,hmax,T)
  h = hmax;
  gt(1) = t0;
  gw(1) = y0;
  t = t0;
  w = y0;

  while (0 == 0)
    k1 = h*funct(t,w);
    k2 = h*funct(t+h/4,w+k1/4);
    k3 = h*funct(t+3*h/8,w+3*k1/32+9*k2/32);
    k4 = h*funct(t+12*h/13,w+1932*k1/2197-7200*k2/2197+7296*k3/2197);
    k5 = h*funct(t+h,w+439*k1/216-8*k2+3680*k3/513-845*k4/4104);
    k6 = h*funct(t+0.5*h,w-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40);

    sigma = abs(k1/360-128*k3/4275-2197*k4/75240+k5/50+2*k6/55);
    if (sigma < T)
      % We accept w as an approximation of y(t) .  w is an approximation
      % of  y(t)  given by a Runge-Kutta method of order four.
      t = t+h;
      w = w+25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
      gt = [gt;t];
      gw = [gw;w];
    end
    
    % We have reached  tf  and the program should stop.
    if (t >= tf)
      return;
    end

    if (sigma == 0)
      % We choose a large value for  q  if the error seems to be negligable.
      q = 5;
    else
      q = (T/sigma)^(0.25);
    end

    % We choose the step size less than haxm and larger than hmin
    % such that the local error should still be less than  T.
    if (q < 0.1)
      % We do not reduce the step-size h to less than 1/10
      % its original size.
      h = 0.1*h;
    elseif (q > 4)
      % We do not increase the step-size h to more than 4 times
      % its original size or hmax.
      h = min(4*h,hmax);
    else
      h = min(h*q,hmax);
    end

    % We make sure than the step-size is not smaller than hmin.
    if (h < hmin)
      break;
      gt = NaN;
      gw = NaN;
    end

    % We adjust the step size if we are going to exceed  tf  at the next
    % step.
    if (t + h > tf)
      h = tf - t;
    end
  end
end
