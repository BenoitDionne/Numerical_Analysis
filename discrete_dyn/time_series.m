%  x = time_series(funct, mu, I, X0, N, UNIT)
%
%  Program to draw the time series (solution) of a discrete
%  dynamical system.
%  It also returns the last value of the iteration.
%
%  Given:
%    The iterative function  funct  that define the discrete
%      dynamical system   x_{n+1} = f(x_n , mu).
%    The parameter  mu  .  it can be a vector if the function accepts
%      a vector as argument.
%    The interval  I  such that  f:I -> I .  the endpoints are
%      are only used to determine the minimal range for the vertical
%      axis
%    The initial condition  X0.
%    The number of iterations  N.
%
function x = time_series(funct, mu, I, X0, N)
   % compute the iterations
   t = 0;
   x = X0;
   X = X0;
   for i=1:N
     x = funct(x,mu);
     t = [ t; i];
     X = [ X; x];
   end

   MIN = min([X;I(1)]);
   MAX = max([X;I(2)]);

   % plot the time series
   clf
   plot(t,X,'b.')
   axis( [ -1, N+1, MIN, MAX ]);
   grid on;
   hold on;
   xlabel('t');
   ylabel('x');
end
