f = @(x) sin(log(x));
i = 0:1:4;
h = 0.8./(2.^i);
col = (f(3+h)-f(3-h))./(2*h);
T = richardson(col);

arrayToLaTeX('table.txt',T,[1:1:9],1,0,'Richarson table','RichTable',["$L_h(f)$"," ","$L_h^1(f)$"," ","$L_h^2(f)$"," ","$L_h^3(f)$"," ","$L_h^4(f)$"], []);

% |L_h^i(f) - L_{2h}^{i-1}(f)|

abs(T(2,3) - T(1,1))
abs(T(3,5) - T(2,3))
abs(T(4,7) - T(3,5))
abs(T(5,9) - T(4,7))
