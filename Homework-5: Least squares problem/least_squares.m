clear all; clc;

format long;

% given
m = 50; % no. of grid points
n = 10; % no. of coefficients
t = linspace(0, 1, m)';  % grid points
b = cos(4*t);  % function values

% Vandermonde matrix A
A = zeros(m, n);
for i = 1:n
    A(:, i) = t.^(i-1); 
end

% Cholesky factorization
R_chol = cholesky(A'*A);
x_chol = forward_backward_sub(R_chol, A'*b);
fprintf('\nCholesky solution: \n');
disp(x_chol);

% Reduced QR factorization using modified Gram-Schmidt
[Q_mgs, R_mgs] = mod_gram_schmidt(A);
x_mgs = back_sub(R_mgs, Q_mgs'*b);
fprintf('\nModified Gram-Schmidt solution: \n');
disp(x_mgs);

% Full QR factorization using Householder triangularization
x_h = householder(A, b);
fprintf('\nHouseholder triangularization solution: \n');
disp(x_h);

% Built-in MATLAB function
x_builtin = A \ b;
fprintf('\nBuilt-in solution: \n');
disp(x_builtin);