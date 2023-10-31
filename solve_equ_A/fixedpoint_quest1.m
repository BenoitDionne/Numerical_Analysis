g = @(x) 5./sqrt(x);

clf
x = 1:0.01:4;
y = g(x);
plot(x,y,'b')
grid on
hold on

plot([1,4],[1,4],'r')
xlabel('x')
ylabel('y')
title('y = 5/x^{1/2}')

K = 5/2^(5/2)
b = -(1/log(K))*log(10^5*abs(3-5/sqrt(2))/(1-K))
n = ceil(b)

y = fixedpoint(g,3,10^(-5),100)
