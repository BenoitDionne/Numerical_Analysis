tol = 10^(-10);

% p(x) = x^3 - 5.974925987 x^2 + 9.734512519 x - 2.617993878
a = [-2.617993878, 9.734512519, -5.974925987, 1];

% Newton with Horner with p
r1 = realroot(a,1,100,tol)

% Deflated polynomial
% q(x) = b_2 x^2 + b_1 x + b_0
b = deflated(a, r1)

% Newton with Horner with q
r2 = realroot(b,1,100,tol)

% Newton with Horner with p
r2 = realroot(a,r2,100,tol)

% Deflated polynomial
% r(x) = b_1 x + b_0
b = deflated(b, r2)

r3 = -b(1)/b(2)

% Newton with Horner with p
r3 = realroot(a,r3,100,tol)
