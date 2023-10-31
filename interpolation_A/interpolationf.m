% y = interpolation(funct, points, x)
%
% Compute the Newton interpolation polynmial of the function
% and evaluate it (using a nested polynomial) at the point x given.
%
% Given:
% - funct : This is the function to be interpolated.  The function must
%           handle vector as arguments.
% - points : Those are the interpolation points.  They must be distinct.
%            A row matrix of values.
% - x : Are the values where the interpolating polynomial will be evaluated.
%       A row of values.
%
% Example;
%   pts = linspace(-1,1,4);
%   f = @(x) abs(x);
%   y = interpolation(f, pts, 0.1);
%
function y = interpolation(funct, points, x)

    % default arguments
    arguments
        funct;
        points (1,:) double;
        x double;
    end

    pts = points';
    X = x';
    N = size(pts,1);
    C = zeros(N,N);
    C(:,1) = funct(pts);
    
    for i=2:1:N
        C(1:1:N-i+1,i) = (C(2:1:N-i+2,i-1) - C(1:1:N-i+1,i-1))./(pts(i:1:N) - pts(1:1:N-i+1));
    end
    
    y = C(1,N);
    for j=1:1:N-1
        y = C(1,N-j) + (X-pts(N-j)).*y;
    end
end
