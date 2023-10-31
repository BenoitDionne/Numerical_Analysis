% [tt,ww] = linearFDM(f,Ai,Aii,Baa,Bab,Bbb,yc,N,a,b)
%
% Solve the boundary value problem y' - A(t) = f(f) with
% B_a y(a) + B_b y(b) = y_c  for $a <= t <= b .
%
% The trapezoidal scheme is used.
%
% Return: The n x (N+1) matrix ww that contains the approximation
%         w_{k,t_i} of y_k(t_i) and the vector tt that contains t_i
%         for 1 \leq k \leq n and 0 <= i <= N.
%
% Given:
%   F: a function from [a,b] x R to R^n,
%   Ci: a function from [a,b] x R to the n x n matrices,
%   Cii: a function from [a,b] x R to the n x n matrices,
%   B_a^{[a]}: a (n-q) x n matrix (Baa in the code below),
%   B_a^{[b]}: a q x n matrix (Bab in the code below),
%   B_b^{[b]}: a q times n matrix (Bbb in the code below),
%   y_c: a vector in R^n (yc in the code below),
%   N > 2: the number of partitions of [a,b], and
%   a and b: the endpoints of the t interval for the domain.
%
function [tt,ww] = linearFDM(F,Ci,Cii,Baa,Bab,Bbb,yc,N,a,b)
    n = length(yc);
    q = size(Bbb,1);
    nmq = n - q;
    h = (b-a)/N;
    
    % We construct the matrix A and the vector F
    A = zeros(n,n,N+1);        % A(:,:,i) = A_{i-1}  for 1 <= i <= N+1
    B = zeros(n,n,N);          % B(:,:,i) = B_i      for 1 <= i <= N
    C = zeros(n,n,N+1);        % C(:,:,i) = C_{i-1}  for 1 <= i <= N+1
    FF = zeros(n,N+1);         % F(:,:,i) = \tilde{f}_{i-1}  for 1 <= i <= N+1

    A(1:nmq,:,1) = Baa;
    C(nmq+1:n,:,N+1) = Bab;
    t = (N-1)*h;
    Civ = Ci(t,h);
    Ciiv = Cii(t+h,h);
    B(1:nmq,:,N) = Civ(q+1:n,:);
    A(1:nmq,:,N+1) = Ciiv(q+1:n,:);
    A(nmq+1:n,:,N+1) = Bbb;
    FF(1:nmq,1) = yc(1:nmq,1);
    tt = [];
    for i=1:1:N
        t = a + (i-1)*h;
        tt = [tt t];
        Civ = Ci(t,h);
        Ciiv = Cii(t+h,h);
        A(nmq+1:n,:,i) = Civ(1:q,:);
        C(nmq+1:n,:,i) = Ciiv(1:q,:);
        B(1:nmq,:,i) = Civ(q+1:n,:);
        A(1:nmq,:,i+1) = Ciiv(q+1:n,:);
        v = F(t,h);
        FF(nmq+1:n,i) = v(1:q,1); 
        FF(1:nmq,i+1) = v(q+1:n,1); 
    end
    FF(nmq+1:n,N+1) = yc(nmq+1:n,1);
    tt = [tt b];

    %%% Testing
    %%% AA = zeros(n*(N+1),n*(N+1));
    %%% FFF = zeros(n*(N+1),1);
    %%% 
    %%% for i=1:1:N+1;
    %%%     icr = (i-1)*n;
    %%%     ir1 = icr + 1;
    %%%     ir2 = icr + n;
    %%%     ic1 = icr - n +1;
    %%%     ic2 = icr;
    %%%     ic3 = icr + 1;
    %%%     ic4 = icr + n;
    %%%     ic5 = ic4 + 1;
    %%%     ic6 = ic4 + n;
    %%%     if ( i > 1)
    %%%         AA(ir1:ir2,ic1:ic2) = B(:,:,i-1);
    %%%     end
    %%%     AA(ir1:ir2,ic3:ic4) = A(:,:,i);
    %%%     if ( i < N+1 )
    %%%         AA(ir1:ir2,ic5:ic6) = C(:,:,i);
    %%%     end
    %%%     FFF(ir1:ir2,1) = FF(:,i);
    %%% end
    %%% AA(ir1:ir2,1:n) = C(:,:,N+1);
    %%% AA;
    %%% FFF;
    %%% WWW = linsolve(AA,FFF);
    %%% WWW([1:2:2*N+1],1)' - exp(tt);
    %%%

    % We construct the matrices L nad U
    Ud = zeros(n,n,N+1);       % Ud(:,:,i) = U_{i-1,i-1} , 1 <= i <= N+1
    Uu = zeros(n,n,N);         % Uu(:,:,i) = U_{i-1,i} , 1 <= i <= N
    Ll = zeros(n,n,N);         % Ll(:,:,i) = L_{i,i-1} , 1 <= i <= N
    Lr = zeros(n,n,N);         % Lr(:,:,i) = L_{N.i-1} ,  1 <= i <= N-1

    Ud(:,:,1) = A(:,:,1);
    Uu(:,:,1) = C(:,:,1);
    % Lr(:,:,1) = C(:,:,N+1)*inv(Ud(:,:,1));
    Lr(:,:,1) = linsolve(Ud(:,:,1)',C(:,:,N+1)')';
    for i=1:1:N-2
        % Ll(:,:,i) = B(:,:,i)*inv(Ud(:,:,i));
        Ll(:,:,i) = linsolve(Ud(:,:,i)',B(:,:,i)')';
        Ud(:,:,i+1) = A(:,:,i+1) - Ll(:,:,i)*Uu(:,:,i);
        Uu(:,:,i+1) = C(:,:,i+1);
        % Lr(:,:,i+1) = -Lr(:,:,i)*Uu(:,:,i)*inv(Ud(:,:,i+1));
        Lr(:,:,i+1) = -linsolve(Ud(:,:,i+1)', Uu(:,:,i)'*Lr(:,:,i)')';
    end
    % Ll(:,:,N-1) = B(:,:,N-1)*inv(Ud(:,:,N-1));
    Ll(:,:,N-1) = linsolve(Ud(:,:,N-1)',B(:,:,N-1)')';
    Ud(:,:,N) = A(:,:,N) - Ll(:,:,N-1)*Uu(:,:,N-1);
    Uu(:,:,N) = C(:,:,N);
    % Ll(:,:,N) = (B(:,:,N) - Lr(:,:,N-1)*Uu(:,:,N-1))*inv(Ud(:,:,N));
    Ll(:,:,N) = linsolve(Ud(:,:,N)',(B(:,:,N) - Lr(:,:,N-1)*Uu(:,:,N-1))')';
    Ud(:,:,N+1) = A(:,:,N+1) - Ll(:,:,N)*Uu(:,:,N);
    
    %%% Testing
    %%% LL = zeros(n*(N+1),n*(N+1));
    %%% irN = N*n;
    %%% irN1 = irN + 1;
    %%% irN2 = irN + n;
    %%% for i=1:1:N+1
    %%%     icr = (i-1)*n;
    %%%     ir1 = icr + 1;
    %%%     ir2 = icr + n;
    %%%     ic1 = icr - n + 1;
    %%%     ic2 = icr;
    %%%     ic3 = icr + 1;
    %%%     ic4 = icr + n;
    %%%     ic5 = ic4 + 1;
    %%%     ic6 = ic4 + n;
    %%%     if ( i > 1 )
    %%%         LL(ir1:ir2,ic1:ic2) = Ll(:,:,i-1);
    %%%     end
    %%%     LL(ir1:ir2,ic3:ic4) = eye(n);
    %%%     UU(ir1:ir2,ic3:ic4) = Ud(:,:,i);
    %%%     if ( i < N+1 )
    %%%         UU(ir1:ir2,ic5:ic6) = Uu(:,:,i);
    %%%     end
    %%%     if ( i < N )
    %%%         LL(irN1:irN2,ic3:ic4) = Lr(:,:,i);
    %%%     end
    %%% end
    %%% LL;
    %%% UU;
    %%% max(max(abs(AA - LL*UU)));
    %%% VV = linsolve(LL,FFF);
    %%% WW = linsolve(UU,VV);
    %%% WW([1:2:2*N+1],1)' - exp(tt);
    %%%
    
    % We now solve the system A W = F
    % First, we solve L V = F
    V = zeros(n,N+1);
    V(:,1) = FF(:,1);
    for i=2:1:N+1
        V(:,i) = FF(:,i) - Ll(:,:,i-1)*V(:,i-1);
    end
    for i=1:1:N-1
        V(:,N+1) = V(:,N+1) - Lr(:,:,i)*V(:,i);
    end
    
    % Second, we solve U W = V
    W = zeros(n,N+1);

    % W(:,N+1) = inv(Ud(:,:,N+1))*V(:,N+1);
    W(:,N+1) = linsolve(Ud(:,:,N+1),V(:,N+1));
    for i=N:-1:1
        % W(:,i) = inv(Ud(:,:,i))*(V(:,i) - Uu(:,:,i)*W(:,i+1));
        W(:,i) = linsolve(Ud(:,:,i),V(:,i) - Uu(:,:,i)*W(:,i+1));
    end
    ww = W;
end
