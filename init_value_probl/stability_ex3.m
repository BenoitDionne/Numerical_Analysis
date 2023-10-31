clf

% Adams-Moulton three-step method
am3s = @(z) -24*(z.^2.*(1-z))./(9*z.^3+19*z.^2-5*z+1);

t = linspace(0,2*pi,360);
xx1 = cos(t) +i*sin(t);
xx2 = 1.2*xx1;
xx3 = 0.8*xx1;

% Absolute stability region for the Adams-Moulton three-step method
z = am3s(xx1);
x = real(z);
y = imag(z);
plot(x,y,'k');
grid on;
hold on;
xlabel('x');
ylabel('y');
% axis equal;

z = am3s(xx2);
x = real(z);
y = imag(z);
plot(x,y,'r');

z = am3s(xx3);
x = real(z);
y = imag(z);
plot(x,y,'b');

% Save the figure
print('stability_ex3','-dpng')
