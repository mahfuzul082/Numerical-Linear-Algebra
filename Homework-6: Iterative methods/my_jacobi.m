function [x, num_iter, time, res] = my_jacobi(A, b, x_0, tol, max_iter)

tic;
n = size(A, 1);
x = x_0;
res = zeros(max_iter, 2);  % column 1: iteration no., column 2: residual

for iter = 1:max_iter
    x_old = x;
    for i = 1:n
        x(i) = ( b(i) - A(i, 1:i-1) * x_old(1:i-1) - A(i, i+1:n) * x_old(i+1:n) ) / A(i,i); % update
    end
    
    % store iteration number and residual
    res(iter, 1) = iter;
    res(iter, 2) = norm(b - A*x);
    
    % convergence check
    if res(iter, 2) <= tol
        fprintf("Prescribed tolerance reached. Final residual = %.2e\n", res(iter, 2));
        time = toc;
        num_iter = iter;
        res = res(1:iter, :);  % only till iterations that satisfy tolerance
        return;
    end
end
% if code reaches here, tolerance not met
time = toc;
fprintf("Prescribed tolerance not reached. Final residual = %.2e\n", res(end, 2));
num_iter = max_iter;
end