f = @(x) (x-2).^(1/4);
a = 3;
b = 5;

col1(1,1) = trapezoid(f,a,b,1);
col1(2,1) = trapezoid(f,a,b,2);
col1(3,1) = trapezoid(f,a,b,4);
col1(4,1) = trapezoid(f,a,b,8);
col1(5,1) = trapezoid(f,a,b,16);

T = richardson(col1);

arrayToLaTeX('table.txt',T,[1:1:9],1,0,'table','RichTable',["$L_h(f)$"," ","$L_h^1(f)$"," ","$L_h^2(f)$"," ","$L_h^3(f)$"," ","$L_h^4(f)$"], []);

% |L_h^i(f) - L_{2h}^{i-1}(f)|

abs(T(2,3) - T(1,1))
abs(T(3,5) - T(2,3))
abs(T(4,7) - T(3,5))
abs(T(5,9) - T(4,7))

