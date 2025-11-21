function x = householder(A, b)

[V, R] = householder_vr(A);
b_tilde = householder_b(V, b);
x = back_sub(R, b_tilde);
end