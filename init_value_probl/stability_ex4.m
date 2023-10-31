clf

% Trapezoidal method
tm = @(z) 2*(z-1)./(z+1);

t = linspace(0,2*pi,360);
xx1 = cos(t) +i*sin(t);
xx2 = 1.2*xx1;
xx3 = 0.8*xx1;

% Absolute stability region for the trapezoidal method
z = tm(xx1);
x = real(z);
y = imag(z);
plot(x,y,'k');
grid on;
hold on;
xlabel('x');
ylabel('y');
% axis equal;

z = tm(xx2);
x = real(z);
y = imag(z);
plot(x,y,'r');

z = tm(xx3);
x = real(z);
y = imag(z);
plot(x,y,'b');

axis([-20 25 -15 15 ]);

% Save the figure
print('stability_ex4','-dpng')
