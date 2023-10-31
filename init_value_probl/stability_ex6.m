clf

% Adams-Bashforth method of order two from the exemple in notes
ab2s = @(z) (z.^2 - 1)./(2*z);

t = linspace(0,2*pi,360);
xx1 = cos(t) +i*sin(t);
xx2 = 1.5*xx1;
xx3 = 0.8*xx1;
xx4 = 1.1*xx1;

% Absolute stability region for the Adams-Bashforth method of order two
z = ab2s(xx1);
x = real(z);
y = imag(z);
plot(x,y,'k');
grid on;
hold on;
xlabel('x');
ylabel('y');
% axis equal;

z = ab2s(xx2);
x = real(z);
y = imag(z);
plot(x,y,'r');

z = ab2s(xx3);
x = real(z);
y = imag(z);
plot(x,y,'b');

z = ab2s(xx4);
x = real(z);
y = imag(z);
plot(x,y,'r');

% Save the figure
print('stability_ex6','-dpng')
