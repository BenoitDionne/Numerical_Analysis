%  C = bezier(P,Qcheck,Qhat,N)
%
%  This program construct the piecewise cubic bezier curve associated to
%    the points P(1,:), P(2,:), ... , P(npoints,:) and the
%    control points Qcheck(i,:) associated to left end point P(i,:) 
%    and Qhat(i+1,:) associated to the right end point P(i+1,:)
%    for i = 1 to npoints-1.
%  The Bezier curve is plotted.
%
%  The following information must be given by the user:
%    The points P(i,:) for i=1 to npoints.
%    The control points Qcheck(i,:) of the left end point P(i,:)
%      for i=1 to npoints-1,
%    The control points Qhat(i,:) of the right end point P(i+1,:)
%      for i=1 to npoints-1,
%    The number N of equaly spaced points of the interval [0,1] to
%      use for plotting the curve.  The default value is 100.
%
%  The following is returned by the function.
%    The graph of the piecewise cubic Bezier curve.
%    The coefficients of the cubic Bezier curve
%      phi_i(t) = a_i t^3 + b_i t^2 + c_i t + d_i
%      used on each interval.
%
function [a,b,c,d] = bezier(P,Qcheck,Qhat,N)
    
  % default arguments
  arguments
    P (:,2) double;
    Qcheck (:,2) double;
    Qhat (:,2) double;
    N int32 = 100;
  end

  npoints = size(P,1);
  figure
  hold on
  grid on

  for i=1:1:npoints-1
    a(i,:) = -P(i,:) + 3*Qcheck(i,:) - 3*Qhat(i,:) + P(i+1,:);
    b(i,:) = 3*P(i,:) - 6*Qcheck(i,:) + 3*Qhat(i,:);
    c(i,:) = -3*P(i,:) + 3*Qcheck(i,:);
    d(i,:) = P(i,:);

    t = linspace(0,1,N)
    x = d(i,1)+(c(i,1)+(b(i,1)+a(i,1).*t).*t).*t;
    y = d(i,2)+(c(i,2)+(b(i,2)+a(i,2).*t).*t).*t;
    plot(x,y,'b');
  end

  xlabel('x')
  ylabel('y')
  title('Piecewise Cubic Bezier Curve')
end


