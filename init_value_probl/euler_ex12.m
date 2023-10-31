cla

% First example
f1 = @(t,y) 0.2*t*y;
[t1,w1] = euler_gr(f1, 1,1.5, 1, 5);

t = 1:0.01:1.5;
y1 = exp(0.1*t.^2 -0.1);

plot(t1,w1,'b');
hold on
plot(t,y1,'k');
grid on
xlabel('t')
ylabel('y')
title('y'' = 0.2 t y');

% Second example
figure
f2 = @(t,y) 2*t*y;
[t2,w2] = euler_gr(f2, 1, 1.5, 1, 5);

t = 1:0.01:1.5;
y2 = exp(t.^2 -1);

plot(t2,w2,'b');
hold on
plot(t,y2,'k');
grid on
xlabel('t')
ylabel('y')
title('y'' = 2 t y');

% euler method
function [ts,ws] = euler_gr(funct, t0, tf, y0, N)
  ts = [t0];
  ws = [y0];
  h = (tf-t0)/N;
  for n=1:1:N
      ws = [ ws , ws(n) + h.*funct(ts(n), ws(n)) ];
      ts = [ ts, t0 + n*h ];
  end
end
