function R = cholesky(A)

% no. of rows of A
n = size(A, 1);

% check symmetry
if ~isequal(A, A')
    error("The input matrix is not symmetric. Cholesky factorization not possible.");
else
    fprintf("The input matrix is symmetric. Checking positive-definiteness...\n");
end

% initialize array for storing determinant of leading principal minors
det_arr = zeros(n, 1);

% check positive-definiteness
for i = 1:n
    det_minor = det(A(1:i, 1:i));
    det_arr(i) = det_minor;
    if det_minor <= 0
        error("The input matrix is symmetric, but not positive definite. " + ...
            "Cholesky factorization not possible.");
    end
end

if all(det_arr > 0)
    fprintf("The input matrix is SPD.\n");
end

fprintf("Determinants of leading principal minors:\n");
disp(det_arr);

% initialize R
R = triu(A);

% modified Gaussian elimination
for k = 1:n
    for i = k+1:n
        m = R(k, i)/R(k, k);
        R(i, i:n) = R(i, i:n) - m*R(k, i:n);
    end
    R(k, k:n) = R(k, k:n)/sqrt(R(k, k));
end
end