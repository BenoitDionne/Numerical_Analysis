d = [[ 0 , 0];[pi/4 , sqrt(2)/2];[pi/4 , sqrt(2)/2];[pi/4 , -sqrt(2)/2];[pi/2 , 1]];
[z, table] = divideddiff(d)

% To write an array to a text file
% writematrix(table,'result.txt','Delimiter','tab')

% To write an array as a LaTeX table to a text file
arrayToLaTeX('result.txt',table,[1,2,3,4,5,6],1,0,'Table of Divided Differences', 'cosdivdiff',["x_i","f[\cdot]","f[\cdot,\cdot]","f[\cdot,\cdot,\cdot]","f[\cdot,\cdot,\cdot,\cdot]","f[\cdot,\cdot,\cdot,\cdot,\cdot]"],[]);

% To write an array to a text filewhile controling the format
% fid = fopen('result.txt','wt');
% if fid > 0
%   fprintf(fid,'%.12f & %.12f & %.12f & %.12f & %.12f & %.12f\n',table');
%   fclose(fid);
% end

clf
x = linspace(0,pi/2,100);
y = polynomial(d(1:4,1)',z,x);
plot(x,y,'b')
grid on
hold on

y = cos(pi/2 - x);
plot(x,y,'r')
xlabel('x')
ylabel('y')

w = polynomial(d(1:4,1)',z,pi/8)
