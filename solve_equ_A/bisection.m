%  x = bisection(funct,a,b,tol)
%   We use the bisection algorithm to find a root
%   of the equation  f(x) = 0.      
%
%   The following must be given:
%     The two endpoints  a  and  b  of the interval  [a,b].
%     The maximal tolerated error  tol.  This value should not
%         be smaller than the smallest number that the computer
%         may handle.
%     The function  f  is given in  funct.
%
%   The program return:
%     NaN  a_0  b_0  NaN             NaN
%     x_1  a_1  b_1  sign of f(x_1)  sign of f(a_0)
%     x_2  a_2  b_2  sign of f(x_2)  sign of f(a_1)
%      ...
%     until  |a_n - b_n| < tol.
%     x = a_n or b_n is then an approximation of a root of f(x) correct to tol.
%
function x = bisection(funct,a,b,tol)
    fa = funct(a);
    fb = funct(b);
    x = NaN;
   
    if ( a >= b )
      disp(['a must be smaller than b.'])
      return;
    end

    % We compute the theoritical number of iterations needed to reach
    % the accuracy requested.  This also prevent any infinite loops.
    %
    % From  (b-a)/(2^n) < tol  we get
    M = ceil(log2((b-a)/tol));

    disp(sprintf('n   x_n   a_n   b_n   sign(f(x_n))   sign(f(a_{n-1}))'));
    data = [0,NaN,a,b,NaN,NaN];
    disp(sprintf('%3d  %0.7f  %0.7f  %0.7f  %2d  %2d',data))

    % We replace fa*fb > 0 by a simple comparison of the signs of these
    % values.  We avoid a multiplication.
    if ( sign(fa) == sign(fb) )
        x = NaN;
          disp(sprintf('The bisection algorithm cannot be used because f(%f) = %f and f(%f) = %f have the same sign.',a,b,fa,fb));
        return;
    end

    p = b - a;
    % We stop at i = N-1 because x_N is computed at i = N-1.
    for i=1:N-1
        p = p/2;

        % Instead of using the formula  (a+b)/2  to compute the middle
        % point, we simply add p to a.
        x = a + p;
        fx = funct(x);

        % The test  fx == 0  is not reliable because it is extremely rare
        % that the numerical evaluation of a function will give exactly 0.
        % We replace this test by  abs(fx) < 2*realmin , where realmin is the
        % smallest number that the computer may handle.
        if ( abs(fx) <= 2*realmin )
            data
            disp(sprintf('After %d iterations, we find %.7f as an approximation of a root of f(x).',i,x));
            return;
        end

        % We replace fa*fx > 0 by a simple comparaison of the signs of these
        % values.  We avoid a multiplication.
        % We also store the value fx of f at the midpoint x into fa if
        % a takes the value x or into fb if b takes the value x.
        % This eliminates the need to compute f again at x.
        sgfx = sign(fx);
        sgfa = sign(fa);
        if ( sgfx ~= sgfa )
            b = x;
            fb = fx;
        else
            a = x;
            fa = fx;
        end

        data = [i,x,a,b,sgfx,sgfa];
        disp(sprintf('%3d  %0.7f  %0.7f  %0.7f  %2d  %2d',data));
    end

    disp(sprintf('After %d iterations, we find %.7f as an approximation of a root of f(x).',i,x));

end

