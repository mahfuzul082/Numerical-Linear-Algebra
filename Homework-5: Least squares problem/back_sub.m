function x = back_sub(R, b_tilde)

%n = length(b_tilde);
n = size(R, 2);
x = zeros(length(b_tilde), 1);

for i = n:-1:1
    if i == n
        x(i) = b_tilde(i) / R(i, i);
    else
        x(i) = (b_tilde(i) - R(i, i+1:n) * x(i+1:n)) / R(i, i);
    end
end
end