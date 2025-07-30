function dxx = CW_eom(xx,~,u,n)

% Build matrices
M1 = zeros(3); M1(1,1) = 3; M1(3,3) = -1;
M2 = zeros(3); M2(2,1) = -2; M2(1,2) = 2;

A = [zeros(3), eye(3); n^2*M1, n*M2];
B = [zeros(3); eye(3)];

% COmpute derivative of state
dxx = A*xx + B*u;