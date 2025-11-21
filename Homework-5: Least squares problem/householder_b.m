function b_tilde = householder_b(V, b)

% no. of rows for b
[m, ~] = size(b);

for k = 1:size(V, 2)
    % extract v_k
    vk = V(k:m, k);
    
    % update b(k:m)
    b(k:m) = b(k:m) - 2 * vk * (vk' * b(k:m));
end

% final updated b
b_tilde = b;
end