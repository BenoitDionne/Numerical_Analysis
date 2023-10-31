
x = -3*pi/2+0.001:0.01:-pi/2-0.001;
y=tan(x);
plot(x,y,'k')
grid on
hold on

x = -pi/2+0.001:0.01:pi/2-0.001;
y=tan(x);
plot(x,y,'k')

x = pi/2+0.001:0.01:3*pi/2-0.001;
y=tan(x);
plot(x,y,'k')

x=-3*pi/2:0.01:3*pi/2;
y = exp(x);
plot(x,y,'b')

% x axis
x = [-3*pi/2 3*pi/2];
y = [0 0];
plot(x,y,'k')

% y axis
x = [0 0];
y = [-20 20];
plot(x,y,'k')

xlabel('x')
ylabel('y')
title('y=tan(x)  and  y = e^x')
axis([-3*pi/2 3*pi/2 -20 20])
