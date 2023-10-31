P = [ [0,2] ; [1,3] ; [3,3] ; [4,2] ; [5,2] ; [5.5, 1] ];
Qcheck = [ [1,2] ; [1.5,3.5] ; [3.5,2.5] ; [4.5,1.5] ; [5.5,1.5] ];
Qhat = [ [0.5,2.5] ; [2.5,3.5] ; [3.5,2.5] ; [4.5,2.5] ; [5,1] ];
N = 100;

[a,b,c,d] = bezier(P,Qcheck,Qhat,N);

for i=1:1:5
  disp(sprintf('\\begin{pmatrix} %f \\\\ %f \\end{pmatrix} t^3 + \\begin{pmatrix} %f \\\\ %f \\end{pmatrix} t^2 + \\begin{pmatrix} %f \\\\ %f \\end{pmatrix} t + \\begin{pmatrix} %f \\\\ %f \\end{pmatrix}', a(i,1), a(i,2), b(i,1), b(i,2), c(i,1), c(i,2), d(i,1), d(i,2)));
end
