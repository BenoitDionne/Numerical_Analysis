% [tt,ww] = shooting(f,A,Ba,Bb,yc,N,a,b)
%
% Solve the boundary value problem y' - A(t) = f(f) with
% B_a y(a) + B_b y(b) = y_c  for $a <= t <= b .
%
% The classical Fourth order Runga-Kutta is use to solve the initial
% value problems in the algorithm.  For N given, the step size is
% h = (b-a)/N  and  t_i = a + i h  for  0 <= i <= N.
%
% Return: The vector ww that contains the approximation w_i of
%         y_i = y(t_i) for 0 <= i <= N and the vector tt that contains
%         t_i for 0 <= i <= N.
%
% Given:
%   f: a function from [a,b] to R^n,
%   A: a function from [a,b] to the n x n matrices,
%   B_a: a n x n matrix (Ba in the code below),
%   B_b: a n times n matrix (Bb in the code below),
%   y_c: a vector in R^n (yc in the code below),
%   N: the number of partitions of [a,b], and
%   a and b: the endpoints of the t interval for the domain.
%
function [tt,ww] = shooting(f,A,Ba,Bb,yc,N,a,b)
    funct1 = @(t,y) A(t)*y + f(t);
    funct2 = @(t,y) A(t)*y; 
    h = (b-a)/N;
    n = length(yc);

    %   Solve the initial value problem
    %   y'(t) - A(t) y(t) = f(t)  with  y(a) = y_c
    [tt,ww1] = rgkt4(funct1,h,N,a,yc);

    %   Solve the initial value problems
    %   y'(t) - A(t) y(t) = 0  with  y(a) = e_i
    %   for 1 <= j <= n
    WW = repmat(NaN,n,N+1,n);
    for j=1:1:n
        yj = zeros(n,1);
        yj(j) = 1;
        [tt,ww2] = rgkt4(funct2,h,N,a,yj);
        WW(:,:,j) = ww2;
    end

    %  Solve yc -B_a y_0(a) - B_b y_0(b) = Q d
    %  with Q = B_a + B_b Y(b)
    Y = yc - Ba*ww1(:,1) -Bb*ww1(:,N+1);
    Q = Ba + Bb*squeeze(WW(:,N+1,:));
    d = linsolve(Q,Y);

    ww2 = repmat(0,n,N+1);
    for j=1:1:n
        ww2 =  ww2 + d(j)*squeeze(WW(:,:,j));
    end
    ww = ww1 + ww2;
end
