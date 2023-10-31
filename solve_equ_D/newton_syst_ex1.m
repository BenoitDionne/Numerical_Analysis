F = @(x) [3*x(1,1).^2-x(2,1).^2 ; 3*x(1,1).*x(2,1).^2-x(1,1).^3-1];
DF = @(x,y) [6*x(1,1) , -2*x(2,1) ; 3*(x(2,1).^2-x(1,1).^2) , 6*x(1,1).*x(2,1) ];
%     Don't leave blank spaces in math formulae used to define functions.
%     Matlab doesn't accept that.

Y = newton_syst(F,DF, [1;1], 100 , 0.000001)
