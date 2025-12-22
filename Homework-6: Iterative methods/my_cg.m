function [x, num_iter, time, res] = my_cg(A, b, tol)

tic;
m = size(A, 1);
max_iter = m;
res = zeros(m, 2);  % column 1: iteration no., column 2: residual

% initial values
x = zeros(m, 1);
r = b;            
p = r;            

r_old = r' * r;

for n = 1:max_iter
    % step length
    Ap = A * p;
    alpha = r_old / (p' * Ap);
    
    % approximate solution
    x = x + alpha * p;
    
    % residual
    r = r - alpha * Ap;
    
    % store iteration number and residual
    res(n, 1) = n;
    res(n, 2) = norm(r);
    
    % convergence check
    if res(n, 2) <= tol
        fprintf("Prescribed residual reached. Final residual = %.2e\n", res(n, 2));
        time = toc;
        num_iter = n;
        res = res(1:n, :);
        return;
    end
    
    % improvement for search direction
    r_new = r' * r;
    beta = r_new / r_old;
    
    % new search direction
    p = r + beta * p;
    
    r_old = r_new;
end
% if code reaches here, tolerance not met
time = toc;
fprintf("Prescribed tolerance not reached. Final residual = %.2e\n", res(end, 2));
num_iter = max_iter;
end