function [Q, R] = mod_gram_schmidt(A)

% size of input A
[m, n] = size(A);

%initializing Q and R
Q = zeros(m, n);
R = zeros(n, n);

% initializing V
V = A;

% loop for computing Q and R
for i = 1:n
    R(i,i) = sqrt(V(:, i)' * V(:, i));
    Q(:, i) = V(:, i) / R(i,i);
    for j = (i+1):n
        R(i,j) = Q(:, i)' * V(:, j);
        V(:, j) = V(:, j) - R(i, j) * Q(:, i);
    end
end
end