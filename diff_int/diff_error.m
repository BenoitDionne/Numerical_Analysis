clf
h = 0.001:0.0001:0.01;
y = 10^(-8)./(2*h) + h.^2/24;
plot(h,y)
grid on
xlabel('h')
ylabel('R')
title('R = 10^{-8}/(2h) + h^2/24')
