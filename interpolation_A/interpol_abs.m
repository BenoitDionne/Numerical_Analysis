clf

f = @(x) abs(x);

% graph of y=f(x) = |x|
aa = [-1,0,1];
bb = [1,0,1];
plot(aa,bb,'b')
hold on
grid on

% Interpolating polynomial at -1 + 2i/3 for i=0,1,2,3
pts = linspace(-1,1,4);
x = -1:0.01:1;
y = interpolation(f,pts,x);
plot(x,y,'k:')

% Interpolating polynomial at -1 + 2i/7 for i=0,1,2,...,7
pts = linspace(-1,1,8);
x = -1:0.01:1;
y = interpolation(f,pts,x);
plot(x,y,'r--')

% Interpolating polynomial at -1 + 2i/15 for i=0,1,2,...,15
pts = linspace(-1,1,16);
x = -1:0.01:1;
y = interpolation(f,pts,x);
plot(x,y,'g-.')

xlabel('x')
ylabel('y')
title('polynomial interpolation of y = |x|')
text(-0.9, 0.5, 'y=p_{15}(x)')
text(0.3, 0.2, 'y=p_7(x)')
text(0.6, 0.45, 'y=p_3(x)')

