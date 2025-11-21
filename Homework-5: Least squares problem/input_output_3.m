clear all; clc;

format shortE;

% input
A = [1 0; 1 1; 1 2; 1 3];

% compute V and R
[V, R] = householder_vr(A);
fprintf("V = \n");
disp(V);
fprintf("R = \n");
disp(R);

% compute Q
Q = householder_q(V);
fprintf("Q = \n");
disp(Q);

% check Q's orthogonality
if norm(Q'*Q - eye(size(Q))) < 1e-15
    fprintf("Q is orthogonal \n");
else
    fprintf("Q is not orthogonal \n");
end

% check whether A can be recovered
A_ = Q*R;
fprintf("\nRecovered A = \n");
disp(A_);