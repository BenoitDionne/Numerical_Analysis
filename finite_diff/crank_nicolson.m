% [X,T,M] = crank_nicolson(f,gb,gl,gr,c,N,M,a,b,t1,t2)
%
% Solve the heat equation with forcing u_t = c^2 u_xx + f on the domain
%   R = [a,b] x [t1,t2]
% with boundary conditions
%   u(a,t) = gl(t) and u(b,t) = gr(t)
% and the initial conditon
%   u(x,t1) = gb(x)
% using
%   deltax = (b-a)/N and deltat = (t2-t1)/M
%
% Return: X contains the x-coordinates x_0, x_1, ..., x_N of the mesh points,
%         T contains the y-coordinates t_0, t_1, ..., t_N of the mesh points,
%         and U contains the approximation of u at the mesh points;
%            U_{i,j} = u(x_{i-1},t_{j-1}) for 1 <= i <= N+1 and
%            1 <= j <= M+1.
%
% Given:
%   f: the right hand side,
%   gb: the value of u for t = t1,
%   gl: the value of u for x = a,
%   gr: the value of u for x = b,
%   N: the number of partitions of [a,b],
%   M: the number of partitions of [t1,t2],
%   a and b: the endpoints of the x interval for the domain, and
%   t1 and t2: the endpoints of the t interval for the domain.
%
function [X,T,U] = crank_nicolson(f,gb,gl,gr,c,N,M,a,b,t1,t2)
  np1 = N + 1;
  mp1 = M + 1;
  U = repmat(NaN,np1,mp1);
  X = linspace(a,b,np1);
  T = linspace(t1,t2,mp1);

  % Initial condition
  for i=1:1:np1
      U(i,1) = gb(X(i));
  end
  % Boundary conditions
  for j=1:1:mp1
      U(1,j) = gl(T(j));
      U(np1,j) = gr(T(j));
  end

  deltax = (b-a)/N;
  deltat = (t2-t1)/M;
  alpha = deltat*(c/deltax)^2/2;

  nm1 = N - 1;
  J = repmat(0,nm1,nm1);
  K = repmat(0,nm1,nm1);
  for i=1:1:nm1
      for j=i-1:1:i+1
          if (i == j)
              K(i,j) = -1 + 2 * alpha;
              J(i,j) = 1 + 2 * alpha;
          elseif ( j > 0 && j < N )
              K(i,j) = -alpha;
              J(i,j) = -alpha;
          end
      end
  end

  % We make use of the fact that the matrix Q is block lower triangular,
  % and only the diagonal and lower diagonal contain non-trivial blocks.
  nm2 = N - 2;
  B = repmat(NaN,nm1,1);
  B(1,1) = U(2,1) + alpha*(U(1,1) -2*U(2,1) + U(3,1) + U(1,2)) ...
           + (f(X(2),T(1)) + f(X(2),T(2)))*deltat/2;
  for k=2:1:nm2
      B(k,1) = U(k+1,1) + alpha*(U(k,1) -2*U(k+1,1) + U(k+2,1)) ...
               + (f(X(k),T(1)) + f(X(k),T(2)))*deltat/2;
  end
  B(nm1,1) = U(N,1) + alpha*(U(nm1,1) -2*U(N,1) + U(np1,1) +U(np1,2)) ...
      + (f(X(N),T(1)) + f(X(N),T(2)))*deltat/2;
  % Solve the system J U_2 = B_1 with matlab library
  U(2:N,2) = linsolve(J,B);
  % Solve the system J U_2 = B_1 with gauss()
  % U(2:N,2) = gauss(J,B,1);
  
  for j=2:1:M
      jp1 = j + 1;
      B(1,1) = alpha*(U(1,j) + U(1,jp1)) ...
               + (f(X(2),T(j)) + f(X(2),T(jp1)))*deltat/2;
      for k=2:1:nm2
          B(k,1) = (f(X(k),T(j)) + f(X(k),T(jp1)))*deltat/2;
      end
      B(nm1,1) = alpha*(U(np1,j) + U(np1,jp1)) ...
          + (f(X(N),T(j)) + f(X(N),T(jp1)))*deltat/2;
      % Solve the system J U_{j+1} = B_j - K U_j with matlab library
      U(2:N,jp1) = linsolve(J,B-K*U(2:N,j));
      % Solve the system J U_{j+1} = B_j - K U_j with gauss()
      % U(2:N,jp1) = gauss(J,B-K*U(2:N,j),1);
  end
end
