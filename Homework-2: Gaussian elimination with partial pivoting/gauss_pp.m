function x = gauss_pp(A, b)

%no. of rows of A
n = size(A, 1);

for k=1:n-1
    %partial Pivoting
    [~, p] = max(abs(A(k:n, k)));
    p = p + k - 1;
    if p ~= k
        %row swapping
        A([k p], :) = A([p k], :);
        b([k p]) = b([p k]);
    end
    %forward elimination
    for i = k+1:n
        m = A(i,k) / A(k,k);
        A(i, k:n) = A(i, k:n) - m * A(k, k:n);
        b(i) = b(i) - m * b(k);
    end
end

x = zeros(n, 1);
%backward substitution
for i = n:-1:1
    x(i) = (b(i) - A(i,i+1:n)*x(i+1:n)) / A(i,i);
end
return
end