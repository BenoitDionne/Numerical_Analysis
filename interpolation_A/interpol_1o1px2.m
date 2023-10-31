clf

f = @(x) 1./(1+x.^2);

% graph of y = f(x) = 1/(1+x^2)
x = -5:0.01:5;
y = f(x);
plot(x,y,'k')
hold on
grid on

% Interpolating polynomial at -5 + i for i=0,1,2,...,10
pts = linspace(-5,5,11);
x = -5:0.01:5;
y = interpolation(f,pts,x);
plot(x,y,'b')

xlabel('x')
ylabel('y')
title('y = 1/(1+x^2) and  p_{10}(x)');
