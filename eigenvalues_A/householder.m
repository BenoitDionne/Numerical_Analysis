%  B = householder(A)
%
%  Householder Reduction Algorithm to produce a matrix B in the Hessemberg
%    form which is conjugate to the given matrix A.
%  Input: The matrix A.
%  Output: The matrix B.
function B = householder(A)
    dim = size(A,1);

    for i = 1:dim-2
        v = zeros(dim,1);
        z = v;

        % s_i
        s = norm(A(i+1:dim,i));

        if ( s ~= 0 )
            % alpha_i
            alpha = 1/( (s + abs(A(i+1,i)))*s );

            % v_i
            v(i+1,1) = sign(A(i+1,i))*s + A(i+1,i);    
            z(i+1,1) = alpha*v(i+1,1);
            v(i+2:dim,1) = A(i+2:dim,i);
            z(i+2:dim,1) = alpha*v(i+2:dim,1);

            % x_i  and  y_i
            x = A*z;
            y = A'*z;

            %  mu_i
            mu = (alpha * (v' *x))/2;

            % p_i  and  q_i
            z = mu*v;
            p = y - z;
            q = x - z;

            %  A_i
            A = A - v * p' - q*v';
        end
    end
    B = A;
end
