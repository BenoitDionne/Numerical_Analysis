% [tt,ww] = par_shooting(f,A,Baa,Bab,Bbb,yc,M,N,a,b)
%
% Solve the boundary value problem y' - A(t) = f(f) with
% B_a y(a) + B_b y(b) = y_c  for $a <= t <= b .
%
% The classical Fourth order Runga-Kutta is use to solve the initial
% value problems in the algorithm.  For M given, the step size is
% h = (t_{i+1}-t_i)/M  and  t_{i,j} = t_i + j h  for  0 <= j <= M.
%
% Return: The n x (N M+1) matrix ww that contains the approximation
%         w_{k,i*M+j} of y_k(t_{i,j}) and the vector tt that contains
%         t_{i*M+j} = t_{i,j} for 1 \leq k \leq n, 0 <= i < N,
%         and 0 <= j < M if i < N-1 or 0 <= j <= M if i = N-1.
%
% Given:
%   f: a function from [a,b] to R^n,
%   A: a function from [a,b] to the n x n matrices,
%   B_a^{[a]}: a (n-q) x n matrix (Baa in the code below),
%   B_a^{[b]}: a q x n matrix (Bab in the code below),
%   B_b^{[b]}: a q times n matrix (Bbb in the code below),
%   y_c: a vector in R^n (yc in the code below),
%   N: the number of partitions of [a,b],
%   M: the number of partitions of [t_i,t_{i+1}], and
%   a and b: the endpoints of the t interval for the domain.
%
function [tt,ww] = par_shooting(f,A,Baa,Bab,Bbb,yc,M,N,a,b)
    funct1 = @(t,y) A(t)*y + f(t);
    funct2 = @(t,y) A(t)*y;
    n = length(yc);
    q = size(Bbb,1);
    nmq = n - q;
    ttt = repmat(NaN,M+1,N);
    WW1 = repmat(NaN,n,M+1,N);
    WWW = repmat(NaN,n,M+1,q,N);
    PPi = repmat(NaN,q,q,N+1);
    FFi = repmat(NaN,n,q,N+1);
    yci = repmat(NaN,n,N+1);
    
    % We choose the matrix D_a such that M_a is invertible
    Da = zeros(q,n);
    s = 1;
    for j=1:1:n
        v = zeros(1,n);
        v(j) = 1;
        MM = [Baa ; v];
        if ( rank(MM) > nmq )
            Da(s,:) = v;
            s = s +1;
        end
        if ( rank(Da) == q )
            break;
        end
    end

    % We find the matrix F_a
    Ma = [Baa ; Da];
    Mainv = inv(Ma);
    Fa = Mainv(:,(q+1):n);

    % We now compute the approximation of y_[i,0}(t_j) and
    % V_i(t_j) for 0 <= i < N and 0 \leq j \leq M, and the
    % F_i and R_i for 0 <= i <= N.
    % Warning: the indices i and j in matlab are shifted by 1
    %          because vectors start with the index 1.
    H = (b-a)/N;
    h = H/M; 
    ti = a;
    FFi(:,:,1) = Fa;
    
    % We solve M_a y_{c,0} = y_c instead of Baa y_{c,0} = y_c^{[a]}
    % to ensure that there is only one solution for Matlab to find. 
    yci(:,1) = linsolve(Ma,yc);
    
    for i=1:1:N+1
        %   Solve the initial value problem
        %   y'(t) - A(t) y(t) = f(t)  with  y_{i,0}(t_i) = y_{c,i}
        if ( i <= N )
            [t,ww1] = rgkt4(funct1,h,M,ti,yci(:,i));
            WW1(:,:,i) = ww1;
        end

        %   Solve the initial value problems
        %   y'(t) - A(t) y(t) = 0  with  y_{i,j}(t_i) = F_i e_j
        %   for 1 <= j <= q
        WW = repmat(NaN,n,M+1,q);
        for j=1:1:q
            yj = zeros(q,1);
            yj(j) = 1;
            y = FFi(:,:,i)*yj;
            [tt,ww2] = rgkt4(funct2,h,M,ti,y);
            WW(:,:,j) = ww2;
        end
        if ( i <= N )
            ttt(:,i) = tt;
            WWW(:,:,:,i) = WW;
        end

        % We choose F_{i+1} and Y_{c,i+1} for the next interval
        Vi = squeeze(WW(:,M+1,:));
        [Fi,R] = par_QR(Vi);
        FFi(:,:,i+1) = Fi;
        PPi(:,:,i) = inv(R);
        if ( i <= N )
            yci(:,i+1) = (eye(n) - Fi*Fi')*ww1(:,M+1);
        end
        ti = ti + H;
    end

    % We now find the vector d_i for 0 <= i < N .
    qN = q*N;
    qNp1 = q*(N+1);
    As = zeros(qNp1,qNp1);
    Bs = zeros(qNp1,1);
    for i=1:1:N
        qi = q*i;
        qim1 = qi - q + 1;
        As(qim1:qi, qim1:qi) = - eye(q);
        As(qim1:qi, qi+1:qi+q) = squeeze(PPi(:,:,i));
        Bs(qim1:qi,1) = PPi(:,:,i)*(FFi(:,:,i+1)')*WW1(:,M+1,i);
    end
    As(qN+1:qNp1, 1:q) = Bab*FFi(:,:,1);
    As(qN+1:qNp1, qN+1:qNp1) = Bbb*FFi(:,:,N+1);
    Bs(qN+1:qNp1,1) = yc(nmq+1:n,1) - Bab*yci(:,1) - Bbb*yci(:,N+1);

    D = linsolve(As,Bs);

    % The results
    ww = [];
    tt = [];
    for i=1:1:N
        w = zeros(n,M+1);
        for j=1:1:q
            w =  w + D((i-1)*q+j)*squeeze(WWW(:,:,j,i));
        end
        if ( i < N )
            ww = [ww, WW1(:,1:M,i) + w(:,1:M)];
            tt = [tt, ttt(1:M,i)'];
        else
            ww = [ww, WW1(:,:,i) + w];
            tt = [tt, ttt(:,i)'];
        end
    end
end
