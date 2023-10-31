f = @(t,y) 1 + y./t;
sol = @(t) t.*log(t);
t0 = 1;
tf = 2;
x0 = 0;

h = sqrt(2)/10^4;
disp(sprintf('h = %0.14f',h));
n = ceil((tf-t0)/h);

result = zeros(n+1,6);
h = (tf-t0)/n;
result(1,1) = 0;
result(1,2) = t0;
result(1,3) = x0; 
result(1,4) = sol(t0);
result(1,5) = 0;
result(1,6) = 0;
for k=1:1:n
    ki = k+1;
    result(ki,1) = k;
    result(ki,2) = t0 + k*h;
    result(ki,3) = result(k,3) + h.*f(result(k,2),result(k,3));
    result(ki,4) = sol(result(ki,2));
    result(ki,5) = result(ki,3) - result(ki,4);
    result(ki,6) = abs(result(ki,5)/result(ki,4));
end

plot(result(:,2),result(:,3),'b');
grid on
hold on
plot(result(:,2),result(:,4),'r');
xlabel('t')
ylabel('y')


disp(sprintf('n = %5d',n));
disp(sprintf('h = %0.14f',h));
disp(sprintf('t = %5d',result(n+1,2)));
disp(sprintf('w_n = %0.7f',result(n+1,3)));
disp(sprintf('y_n = %0.7f',result(n+1,4)));
disp(sprintf('abs. error = %0.7f',result(n+1,5)));
disp(sprintf('rel. error = %0.7f',result(n+1,6)));
