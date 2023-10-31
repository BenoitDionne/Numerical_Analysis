tol = 10^(-9);

% p(x) = x^4 - 2 x^3 - 12 x^2 +16 x -40
a = [-40, 16, -12, -2, 1 ];

% Newton with Horner with p
r1 = realroot(a,1,100,tol)

% Deflated polynomial
% q(x) = b_3 x^3 + b_2 x^2 + b_1 x + b_0
b = deflated(a, r1)

% Newton with Horner with q
r2 = realroot(b,1,100,tol)

% Newton with Horner with p
r2 = realroot(a,r2,100,tol)

% Deflated polynomial
% q(x) = c_2 x^2 + c_1 x + c_0
c = deflated(b, r2)

% The last roots are complex conjugated
d = c(2)^2 - 4*c(1)*c(3)

r3 = (-c(2) + sqrt(abs(d))*i)/2
r4 = (-c(2) - sqrt(abs(d))*i)/2

% Newton with Horner with p
r3 = realroot(a,r3,100,tol)
r4 = realroot(a,r4,100,tol)

% Check with Matlab
roots(a(5:-1:1))
