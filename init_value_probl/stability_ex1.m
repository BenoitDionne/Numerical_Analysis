clf

ff = @(t,x) 1 - 2*x;
ss = @(t) (1+exp(-2*t))/2;
t0=0; t1=6; N=300; x0 = 1;

% Euler
result_eul = euler_ex(ff, ss, t0, t1, N, x0);

% Adams Bashforth two steps
result_ada = adams_bashforth_2s(ff, ss, t0, t1, N, x0);

t=linspace(t0,t1,N+1);

subplot(2,2,1);
plot(t,result_eul(:,3)','k');
grid on
xlabel('t')
ylabel('y')
title('Euler Method')

subplot(2,2,2);
plot(t,result_ada(:,3)','k');
grid on
xlabel('t')
ylabel('y')
title('Adams Bashforth 2-Step')

subplot(2,2,3);
plot(t,result_ada(:,4)','k');
grid on
xlabel('t')
ylabel('y')
title('Exact solution')

% Save the figure
print('stability_ex1','-dpng')
