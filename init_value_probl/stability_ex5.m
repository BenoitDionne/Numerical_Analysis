clf

% Runge-Kutta method of order 2
rko2p = @(z) -1+sqrt(2*z -1);
rko2m = @(z) -1-sqrt(2*z -1);

% We must go from -pi to pi instead of 0 to 2pi to avoid the jump of abs()
% at -1 that create a close curve.  abs() in matlab return only one of the
% square roots.  
t = linspace(-pi,pi,360);
xx1 = cos(t)+i*sin(t);
xx2 = 1.2*xx1;
xx3 = 0.8*xx1;

% Absolute stability region for the Adams-Moulton three-step method
z = rko2p(xx1);
x = real(z);
y = imag(z);
plot(x,y,'k');
grid on;
hold on;
xlabel('x');
ylabel('y');
% axis equal;

z = rko2p(xx2);
x = real(z);
y = imag(z);
plot(x,y,'r');

z = rko2p(xx3);
x = real(z);
y = imag(z);
plot(x,y,'b');

z = rko2m(xx1);
x = real(z);
y = imag(z);
plot(x,y,'k');

z = rko2m(xx2);
x = real(z);
y = imag(z);
plot(x,y,'r');

z = rko2m(xx3);
x = real(z);
y = imag(z);
plot(x,y,'b');

% Save the figure
print('stability_ex5','-dpng')

% Another approach

% contour plot for |1 + z + z^2/2|
% ff = @(x,y) sqrt( (1 + x + (x.^2-y.^2)/2).^2 + (y.*(1+x)).^2); 

% x = -2.5:0.01:0.5;
% y = -2:0.01:2;
% [X,Y] = meshgrid(x,y);
% Z = ff(X,Y);
% [M,c] = contour(X,Y,Z,[0.8 1.2]);
% grid on;
% hold on;
% c.ShowText = 'on';
% c.LineColor = 'blue';

% [M,c] = contour(X,Y,Z,[1 1]);
% c.ShowText = 'on';
% c.LineColor = 'black';

% xlabel('x');
% ylabel('y');
% axis equal;
