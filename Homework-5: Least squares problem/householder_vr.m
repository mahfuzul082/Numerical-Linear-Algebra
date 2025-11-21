function [V, R] = householder_vr(A)

% size of A
[m, n] = size(A);

% initialize V
V = zeros(m, n);

for k = 1:n
    % extract the k-th column
    x = A(k:m, k);
    
    % define e1
    e1 = zeros(size(x));
    e1(1) = 1;
    
    % sign of the first component of x
    if x(1) == 0
        s = 1;
    else
        s = sign(x(1));
    end
    
    % compute v_k
    vk = s * norm(x) * e1 + x;
    vk = vk / norm(vk); 
    
    % store vk in V
    V(k:m, k) = vk;
    
    % update of A
    A(k:m, k:n) = A(k:m, k:n) - 2 * vk * (vk' * A(k:m, k:n));
end

% final upper triangular matrix
R = A;
end