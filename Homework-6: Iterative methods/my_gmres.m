function [x, num_iter, time, res] = my_gmres(A, b, tol)

tic;
m = size(A, 1);
max_iter = m;
res = zeros(max_iter, 2);  % column 1: iteration no., column 2: residual

% initial residual
r0_norm = norm(b);

if r0_norm == 0
    fprintf("Initial guess is exact solution.\n");
    x = zeros(m, 1);
    num_iter = 0;
    time = toc;
    res = [0, 0];
    return;
end

Q = zeros(m, max_iter+1);
Q(:,1) = b / r0_norm; % q_1

% Hessenberg matrix
H = zeros(max_iter+1, max_iter);

for n = 1:max_iter
    % Arnoldi iteration
    v = A * Q(:,n);
    for j = 1:n
        H(j,n) = Q(:,j)' * v;
        v = v - H(j,n) * Q(:,j);
    end
    H(n+1,n) = norm(v);
    
    if H(n+1,n) ~= 0
        Q(:,n+1) = v / H(n+1,n);
    end
    
    % least-squares solve using Householder QR
    H_tilde = H(1:n+1, 1:n);
    [V_h, R] = householder_vr(H_tilde);
    Q_hat = householder_q(V_h);
    
    % solve the upper triangular system
    rhs = zeros(n+1, 1);
    rhs(1) = r0_norm;
    rhs_hat = Q_hat' * rhs;
    y = back_sub(R(1:n, 1:n), rhs_hat(1:n));
    
    x = Q(:,1:n) * y; % update
    
    % store iteration number and residual
    res(n, 1) = n;
    res(n, 2) = norm(b - A*x);
    
    % convergence check
    if res(n, 2) <= tol
        fprintf("Prescribed tolerance reached. Final residual = %.2e\n", res(n, 2));
        time = toc;
        num_iter = n;
        res = res(1:n, :);
        return;
    end
end
% if code reaches here, tolerance not met
time = toc;
fprintf("Prescribed tolerance not reached. Final residual = %.2e\n", res(end, 2));
num_iter = max_iter;
end