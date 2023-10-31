%  sum = simpson_adapt(funct,a,b,T,Max)
%
%  We use an adaptive method with the composite Simpson's rule to
%  approximate the integral of a function.
%
%  The following must be given:
%    The function given in  "funct",
%    The endpoints  a  and  b  of the interval  [a,b],
%    The total accuracy  T  required.
%    The maximal number Max of time that an interval can be subdivided.
%
%  The program gives an approximation  sum  of the integral.
%
function sum = simpson_adapt(funct,a,b,T,Max)
  sum = nested_adaptive(funct, a, b, T, Max, simpsonNC(funct,a,b));
end

function sum = nested_adaptive(funct,a,b,T,Max, S)
  if (Max < 0 )
    sum = NaN;
    return;
  end

  mid = (a+b)/2;
  sum1 = S;
  sum2L = simpsonNC(funct,a,mid);
  sum2R = simpsonNC(funct,mid,b);
  sum2 = sum2L + sum2R;
 
  if ( abs(sum1-sum2)/15  < T)
    sum = sum2;
  else
    sum = nested_adaptive(funct,a,mid,T/2,Max-1,sum2L) + ...
          nested_adaptive(funct,mid,b,T/2,Max-1,sum2R);
  end
end

function sum = simpsonNC(funct,a,b)
  sum = (b-a).*(funct(a) + 4*funct((a+b)/2) + funct(b))/6;
end
