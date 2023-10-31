f = @(x) [ x(1,1)^3+(x(1,1)^2)*x(2,1)-x(1,1)*x(3,1)+6 ; exp(x(1,1))+exp(x(2,1))-x(3,1) ; x(2,1)^2-2*x(1,1)*x(3,1)-4 ];
df = @(x) [ [3*x(1,1)^2+2*x(1,1)*x(2,1)-x(3,1) , x(1,1)^2 , -x(1,1)] ; [exp(x(1,1)) , exp(x(2,1)) , -1] ; [-2*x(3,1) , 2*x(2,1) , -2*x(1,1) ]];
x = [ 1 ; 1 ; 1];

y = newton_syst(f,df,x,100,0.000001)