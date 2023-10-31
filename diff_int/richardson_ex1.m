f = @(x) sin(x);
i = 0:1:7;
h = 1.6./(2.^i);
col = (f(1+h)-f(1-h))./(2*h);
T = richardson(col);

arrayToLaTeX('table.txt',T,[1:1:15],1,0,'Richarson table to approximate $\sin(x)$ near $x=1$','RichTable',["$0$"," ","$1$"," ","$2$"," ","$3$"," ","$4$"," ","$5$"," ","$6$"," ","$7$"], []);
