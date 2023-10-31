%  [tt,ww] = rgkt4(funct,h,t,w)
%
%  We use Runge-Kutta method of order four to generate approximations
%  w(i)  of  y( t(i) ),  where  t(i) = t(1) + (i-1)*h 
%  for  i = 2, 3, ...
%
%  The following must be given:
%    The initial conditions  t(1) = t0  and  w(1) = y0,
%    The step size  h.
%    The number of steps N.
%    The function  f(t,y) must be given in the function funct.
%    
%  The program gives:
%     t(i) , the apxroximation  w(i)  of  y( t(i) ) ,
%
function [tt,ww] = rgkt4(funct,h,N,t0,y0)
  tt(1) = t0;
  ww(:,1) = y0;
  h2 = h/2;
  for j=1:N
    tt(j+1) = tt(1)+j*h;
    k1 = h*funct(tt(j),ww(:,j));
    k2 = h*funct(tt(j)+h2,ww(:,j)+k1/2);
    k3 = h*funct(tt(j)+h2,ww(:,j)+k2/2);
    k4 = h*funct(tt(j+1),ww(:,j)+k3);
    ww(:,j+1) = ww(:,j) + (k1+2*(k2+k3)+k4)/6;
  end
end
