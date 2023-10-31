% circle(t0,n)
%   Draw the points (cos(j*t0),sin(j*t0)) for 0 <= j <= n
%   Use to illustrate the fact the those points are dense on
%   the uniit circle when n goes to infinity and t0/pi in not
%   rational.
%
function circle(t0,n)
    t = 0:1:n;
    d = 2*pi;
    r = rem(t*t0,d);
    x = cos(r);
    y = sin(r);
    
    clf;
    plot(x,y,'*')
    grid on
    hold on
    xlabel('x')
    ylabel('y')
end
