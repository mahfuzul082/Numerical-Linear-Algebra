function Q = householder_q(V)

% size of V
[m, n] = size(V);

% initialize Q
Q = eye(m);

for k = n:-1:1
    % extract v_k
    vk = V(k:m, k);

    % transform to Q
    Q(k:m, k:m) = Q(k:m, k:m) - 2 * vk * (vk' * Q(k:m, k:m));
end
end