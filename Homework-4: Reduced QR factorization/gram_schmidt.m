function [Q, R] = gram_schmidt(A)

% size of input A
[m, n] = size(A);

%initializing Q and R
Q = zeros(m, n);
R = zeros(n, n);

% initializing V
V = A;

% loop for computing Q and R
for j = 1:n
    for i = 1:(j-1)
        R(i,j) = Q(:, i)' * A(:, j);
        V(:, j) = V(:, j) - R(i, j) * Q(:, i);
    end
    R(j,j) = sqrt(V(:, j)' * V(:, j));
    Q(:, j) = V(:, j) / R(j,j);
end
end