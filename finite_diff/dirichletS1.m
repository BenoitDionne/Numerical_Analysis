% [X,Y,M] = dirichletS1(f,gb,gt,gl,gr,N,M,a,b,c,d)
%
% Solve the Dirichlet equation  u_xx + u_yy = f on the domain
%   R = [a,b] x [c,d]
% with boundary conditions
%   u(x,c) = gb(x), u(x,d) = gt(x), u(a,y) = gl(y) and u(b,y) = gr(y)
% using
%   deltax = (b-a)/N and deltay = (d-c)/M
%
% Return: X contains the x-coordinates x_0, x_1, ..., x_N of the mesh points,
%         Y contains the y-coordinates y_0, y_1, ..., y_N of the mesh points,
%         and U contains the approximation of u at the mesh points;
%            U_{i,j} = u(x_{i-1},y_{j-1}) for 1 <= i <= N+1 and
%            1 <= j <= M+1.
%
% Given:
%   f: the right hand side,
%   gb: the value of u for y = c,
%   gt: the value of u for y = d,
%   gl: the value of u for x = a,
%   gr: the value of u for x = b,
%   N: the number of partitions of [a,b],
%   M: the number of partitions of [c,d],
%   a and b: the endpoints of the x interval for the domain, and
%   c and d: the endpoints of the y interval for the domain.
%
function [X,Y,U] = dirichletS1(f,gb,gt,gl,gr,N,M,a,b,c,d)
  np1 = N + 1;
  mp1 = M + 1;
  U = repmat(NaN,np1,mp1);
  X = linspace(a,b,np1);
  Y = linspace(c,d,mp1);

  % Boundary conditions
  for i=1:1:np1
      U(i,1) = gb(X(i));
      U(i,mp1) = gt(X(i));
  end
  for j=1:1:mp1
      U(1,j) = gl(Y(j));
      U(np1,j) = gr(Y(j));
  end

  deltax = (b-a)/N;
  deltay = (d-c)/M;
  alpha = (deltay/deltax)^2;

  nm1 = N - 1;
  J = repmat(0,nm1,nm1);
  for i=1:1:nm1
      for j=i-1:1:i+1
          if (i == j)
              J(i,j) = -2 -2 * alpha; 
          elseif ( j > 0 && j < N )
              J(i,j) = alpha;
          end
      end
  end

  mm1 = M - 1;
  nm2 = N - 2;
  MNm1 = mm1*nm1;
  deltay2 = deltay^2;
  Q = repmat(0,MNm1,MNm1);
  B = repmat(NaN,MNm1,1);
  for j=1:1:mm1
      for i=j-1:1:j+1
          if (i == j)
              Q((i-1)*nm1+1:i*nm1,(j-1)*nm1+1:j*nm1) = J; 
          elseif ( i > 0 && i < M )
              Q((i-1)*nm1+1:i*nm1,(j-1)*nm1+1:j*nm1) = eye(nm1);
          end
      end
      jm1 = j - 1;
      jp1 = j + 1;
      if (j == 1)
          B(jm1*nm1+1,1) = -U(2,1) -alpha*U(1,jp1) +f(X(2),Y(jp1))*deltay2;
          for k=2:1:nm2
              B(jm1*nm1+k,1) = -U(k+1,1) + f(X(k+1),Y(jp1))*deltay2;
          end
          B(j*nm1,1) = -U(N,1) -alpha*U(np1,jp1) +f(X(N),Y(jp1))*deltay2;
      elseif (j == mm1)
          B(jm1*nm1+1,1) = -U(2,mp1) -alpha*U(1,jp1) +f(X(2),Y(jp1))*deltay2;
          for k=2:1:nm2
              B(jm1*nm1+k,1) = -U(k+1,mp1) + f(X(k+1),Y(jp1))*deltay2;
          end
          B(j*nm1,1) = -U(N,mp1) -alpha*U(np1,jp1) +f(X(N),Y(jp1))*deltay2;
      else
          B(jm1*nm1+1,1) = -alpha*U(1,jp1) +f(X(2),Y(jp1))*deltay2;
          for k=2:1:nm2
              B(jm1*nm1+k,1) = f(X(k+1),Y(jp1))*deltay2;
          end
          B(j*nm1,1) = -alpha*U(np1,jp1) +f(X(N),Y(jp1))*deltay2;
      end
  end

  % Solve the system Q UU = B with matlab library
  %  UU = linsolve(Q,B);

  % Solve the system Q UU = B with gauss()
  UU = gauss(Q,B,1);

  % Transfer UU to U
  for j=1:1:mm1
      jm1 = j - 1;
      jp1 = j + 1;
      U(2:N,jp1) = UU(jm1*nm1+1:j*nm1,1);
  end
end
