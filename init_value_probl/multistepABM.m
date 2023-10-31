%  [gt,wt] = multistepABM(funct,t0,y0,tf,hmin,hmax,T)
%
%  This program uses Adams-Moulton predictor-corrector method
%  with variable step size to solve an ODE
%    y' = f(t,y)
%  with initial condition y(t0) = y0.
%
%  The following must be given:
%    The maximal step size hmax,
%    The minimal step size hmin,
%    The maximal tolerated error T,
%    The initial conditions y0 at t0,
%    The final value  tf  of  t  where to stop the integration.
%    The function  f(t,y)  must be given in  funct.
%
%  There are two flags: rkflag = 1 if the previous stage used
%    Runge-Kutta method and rkflag = 0  otherwise,
%    last = 1 if we have reached tf or h < hmin at some point
%    and last = 0  otherwise.
%
%  The program gives:
%    t(i) , the apxroximation  w(i)  of  y( t(i) ) ,
%
function [gt,gw] = multistepABM(funct,t0,y0,tf,hmin,hmax,T)
  % last = 1  if we have reached  tf  or  h < hmin  at some point,
  % and  last = 0  otherwise.
  last = 0;
  h = hmin;
  gt(1) = t0;
  gw(1) = y0;
  t(1) = t0;
  w(1) = y0;

  % Given t0 and y0, we use Runge-Kutta of order 4, to compute
  % an approximation w(i+1) of y(t0+i*h) for i = 1, 2 and 3 .
  % The code for the function rgkt4() is given separately.
  [t,w] = rgkt4(funct,h,3,t(1),w(1));

  % rkflag = 1  if the last stage used Runge-Kutta of order 4 and
  % rkflag = 0  otherwise,
  rkflag = 1;
  i=1:4;
  f(i) = funct(t(i),w(i));
  
  while (1==1)
    t(5) = t(4) + h;
    % We use the predictor-corrector method
    predict = w(4) + h*(55*f(4) - 59*f(3) + 37*f(2) - 9*f(1))/24;
    f(5) = funct(t(5),predict);
    correct = w(4) + h*(9*f(5) + 19*f(4) - 5*f(3) + f(2))/24;
    sigma = 19*abs(predict-correct)/(270*h);
    if (sigma < T)
      w(5) = correct;
      f(5) = funct(t(5),correct);
      j=1:4;
      t(j) = t(j+1);
      w(j) = w(j+1);
      f(j) = f(j+1);

      if (rkflag==1)
        %  We accept the three values obtained from Runge-Kutta of order 4
        %  and the one obtained with the predictor-corrector method.
        for j=1:4
          gt = [gt,t(j)];
          gw = [gw,w(j)];
        end
      else
        % We accept the new value obtained with the predictor-corrector
        % method.  The other values have already been accepted.
        % It is at least the second time in a row that we apply
        % the predictor-corrector method.
        gt = [gt,t(4)];
        gw = [gw,w(4)];
      end

      if (last == 1)
        break;
      end
      
      % We have now executed at least one iteration of the
      % predictor-corrector method
      rkflag = 0;
      
      if ( (t(4)+h > tf) || (sigma < T/2) )
        % We now choose a bigger step size.
        if (sigma == 0)
          q = 4.0;
        else
          q = (T/sigma)^0.25;
        end
        h = min([hmax,4*h,q*h]);

        % We check that after the next stage  t  will not exceed  tf.
        if (t(4) + h > tf)
          %  We divide by four because we are now going to use
          %  fourstepsrk and one step of Adams-Moulton to
          %  complete the integration;  we must therefore
          %  have that t(4) + 4*h = tf. 
          h = (tf-t(4))/4;
          last = 1;
        end

        t(1) = t(4);
        w(1) = w(4);
        f(1) = f(4);
        [t,w] = rgkt4(funct,h,3,t(1),w(1));
        rkflag = 1;
        j=2:4;
        f(j) = funct(t(j),w(j));
      end
    else
      % We choose a smaller step size.
      q = max([0.1, (T/sigma)^0.25]);
      h = q*h;
      
      if (h < hmin)
        gt = NaN;
        gw = NaN;
        break
      end

      % We start Runge-Kutta with  t(4)  and  w(4)  if
      % we have used the predictor-corrector method at the previous
      % stage.
      if (rkflag == 0)
        t(1) = t(4);
        w(1) = w(4);
        f(1) = f(4);
      end
      [t,w] = rgkt4(funct,h,3,t(1),w(1));
      rkflag = 1;
      last = 0;
      j=2:4;
      f(j) = feval(funct,t(j),w(j));
    end
  end  
end
