clf

% graph of (x-x_0)(x-x_1)...(x-x_7) for x_i = -5 - 10i/7
xi = linspace(-5,5,8);
x = -5:0.01:5;
y = 1;
for i = 1:8
  y = y.*(x-xi(i));
end
plot(x,abs(y),'b');
hold on
grid on

% graph of (x-x_0)(x-x_1)...(x-x_7) for x_i = -5 cos((2i+1)pi/16)
i = 0:1:7;
xi = 5*cos((2*i+1)*pi/16);
y = 1;
for i = 1:8
  y = y.*(x-xi(i));
end
plot(x,abs(y),'k');
xlabel('x')
ylabel('y')
