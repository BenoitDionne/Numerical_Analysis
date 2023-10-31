cla;
axis equal;
grid on;
hold on;
xlabel('\mu');
ylabel('x');
title('Period doubling route to chaos');

a = 2.8:0.001:4;
s = size(a,2);
x = 0.5*ones(1,s);

for i=1:1:1000
    x = a.*x.*(1-x);
end

for i=1:1:1000
    x = a.*x.*(1-x);
    plot(a,x,'b.','markersize',1);
end
