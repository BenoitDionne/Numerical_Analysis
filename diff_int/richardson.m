% T = richardson(col1)
%
% Generate the Richardson table associated to the date given in the
% first column.
%
% Given:
%   col1 : The first column (a column vector) of the table.
%
function T = richardson(col1)
  % default arguments
  arguments
    col1 (:,1) double;
  end

  N = size(col1,1);
  T = repmat(NaN,N,2*N-1);
  T(:,1) = col1;
    
  for i=2:1:N
    T(i:1:N,2*i-1) = ((4^(i-1))*T(i:1:N,2*i-3) - T(i-1:1:N-1,2*i-3))/(4^(i-1)-1);
    if ( i < N )
      T(i:1:N-1,2*(i-1)) = (T(i:1:N-1,2*i-3) - T(i-1:1:N-2,2*i-3))./(T(i+1:1:N,2*i-3)-T(i:1:N-1,2*i-3));
    end
        %  To get the table aligned on the top
        %
        %        T(1:1:N-i+1,2*i-1) = ((4^(i-1))*T(2:1:N-i+2,2*i-3) - T(1:1:N-i+1,2*i-3))/(4^(i-1)-1);
        %       if ( i < N )
        %           T(2:1:N-i+1,2*(i-1)) = (T(2:1:N-i+1,2*i-3) - T(1:1:N-i,2*i-3))./(T(3:1:N-i+2,2*i-3)-T(2:1:N-i+1,2*i-3));
        %       end
  end
end

