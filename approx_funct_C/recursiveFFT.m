% function Z = recursiveFFT(X)
%
% To compute the Fast Fourier Transform of a N-periodic function x,
% where N = 2^m.
%
% The following must be given
%   The column vector X representing x(j) for j=0,1,2,...,N-1.
%   Because of the periodicity of the function x, only these values
%   are needed.
%
% The function gives
%   A vector Z with the values of the Fast Fourier Transform of x
%   at 0,1,2,...,N-1.  As for the input, only these values are given
%   because of the Fast Fourier Transform of x is N-periodic.
%
function Z = recursiveFFT(X)
    N = length(X);
    if N == 1
        Z = X;
    else
        % We compute the Fourier Transform for x_{2k}
        Y1 = recursiveFFT( X(1:2:N) );

        % We compute the Fast Fourier Transform for x_{1+2k}
        Y2 = recursiveFFT( X(2:2:N) );

        a = [0:N/2-1]';
        Y3 = Y2.*exp(-(2*pi*i)*a/N);
        Z = [Y1+Y3 ; Y1-Y3];
    end
end
