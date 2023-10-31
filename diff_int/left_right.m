%  [l,r] = left_right(funct,a,b,N)
%
%  Compute the Riemann left sum
%    sum_{i=0}^{i=N-1} f(a + i h) h
%
%  and the Riemann right sum
%    sum_{i=1}^{i=N} f(a + i h) h
%
%  where  h = (b-a)/N
%
function [sl,sr] = left_right(funct,a,b,N)
    h = (b-a)/N;
    x = linspace(a,b,N+1);
    r = x(2:N+1);
    l = x(1:N);
    sr = sum(funct(r))*h;
    sl = sum(funct(l))*h;
end
