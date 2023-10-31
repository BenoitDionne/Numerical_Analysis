% function y = Bspline(k,i,t,T)
%
% Function to draw a B-spline
%
% The following information must be given:
%   The degree k of the B-spline
%   The index i of the B-spline
%   The value t where to evaluate the B-spline
%   The knots  T(i) for 1 <= i <= i + k + 1 at least to get the
%      non-null part of the graph of the B-spline
%
% The function return the value of the B-spline at t
%
function y = Bspline(k,i,t,T)
    if ( k > 0 )
        y = (t-T(i))./(T(i+k)-T(i)).*Bspline(k-1,i,t,T) + ...
            (1 - (t-T(i+1))./(T(i+k+1)-T(i+1)) ).*Bspline(k-1,i+1,t,T);
    else
        y = ( t >= T(i) & t < T(i+1) );
    end
end
