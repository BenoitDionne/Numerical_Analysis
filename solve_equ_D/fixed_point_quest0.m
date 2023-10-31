g = @(x) [ (x(2,1) - x(1,1).^3 + 3)/12 ; (2*x(1,1) +x(2,1).^3 +2)/12 ];
x = [ 0 ; 0 ];

fixed_point_syst(g,x,100,0.00001,2)

t2 = -0.2:0.01:0.7;
x2 = t2.^3 + 12*t2 -3;
t1 = -4:0.01:4;
x1 = (-t1.^3 +12*t1-2)./2;

clf
plot(t2,x2,"b")
hold on
grid on
plot(x1,t1,"r")
xlabel('x_1')
ylabel('x_2')
