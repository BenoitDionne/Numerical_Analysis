clf

% Adams-Bashforth four-step method
ab4s = @(z) -24*(z.^3.*(1-z))./(55*z.^3-59*z.^2+37*z-9);

t = linspace(0,2*pi,360);
xx1 = cos(t) +i*sin(t);
xx2 = 1.2*xx1;
xx3 = 0.8*xx1;

% Absolute stability region for the Adams-Bashforth four-step method
z = ab4s(xx1);
x = real(z);
y = imag(z);
plot(x,y,'k');
grid on;
hold on;
xlabel('x');
ylabel('y');
% axis equal;

z = ab4s(xx2);
x = real(z);
y = imag(z);
plot(x,y,'r');

z = ab4s(xx3);
x = real(z);
y = imag(z);
plot(x,y,'b');

% Save the figure
print('stability_ex2','-dpng')
