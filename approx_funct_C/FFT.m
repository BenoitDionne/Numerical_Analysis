% function Z = FFT(X,MP)
%
% To compute the Discrete Fourier Transform of a N-periodic function x.
% We assume that N = P1*P2*...*Pm.
%
% The following must be given
%   The vector MP = [P1,P2, ...Pm] and the column vector X
%   representing x(j) for j=0,1,2,...,N-1.  Because of the periodicity
%   of the function x, only these values are needed.
%
% The function gives
%   A vector Z with the values of the Fast Fourier Transform of x
%   at 0,1,2,...,N-1.  As for the input, only these values are given
%   because of the Fast Fourier Transform of x is N-periodic.
%
function Z = FFT(X,MP)
    m = length(MP);

    % For k = m, Y = X
    Y(:,1) = X;
    for k = m:-1:1
        if k < m
            B = prod(MP(k+1:m));
        else 
            B = 1;
        end
        P = MP(k);
        if ( k > 1 )
            A = prod(MP(1:k-1));
        else
            A = 1;
        end

        x = 1;
        omega = exp(-(2*pi*i)/(P*B))
        p = [0:P-1]';
        for b = 0:B-1
            omega_p =  omega.^(b+B*p);
            for a = 0:A-1

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Do not forget that the column vectors in Matlab are
                % indexed (1,1),(2,1), ...
                Ytempo = ones(P,1)*Y(1+b+B*a+A*B*(P-1),1);
                for s = P-2:-1:0
                    Ytempo = Ytempo.*omega_p + Y(1+b+B*a+A*B*s,1);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % We transfer the information to Z for the next value of m.
                for s = 0:P-1
                    Z(1+b+B*s+B*P*a,1) = Ytempo(s+1);
                end
            end
        end

        % The value of Y for the next value of m.
        Y = Z;
    end

    % The final result
    % For m = 1, Z is the Discrete Fourier Transform of X. 
end
