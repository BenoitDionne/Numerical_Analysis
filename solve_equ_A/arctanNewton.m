cla

x=-6:0.1:6;
y=atan(x);
plot(x,y,'k')
hold on
grid on

z=2*x.*(1+x.^2).^(-1);
plot(x,z,'b')

% x axis
x = [-6 6];
y = [0 0];
plot(x,y,'k')

% y axis
x = [0 0];
y = [-1.5 1.5];
plot(x,y,'k')

xlabel('x')
ylabel('y')
title('y=arctan(x)  and  y=2x/(1+x^2)')
