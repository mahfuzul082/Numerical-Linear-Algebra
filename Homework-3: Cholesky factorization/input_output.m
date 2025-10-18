clear all; clc;

% inputs
A = [4 1 1 1; 1 3 -1 1; 1 -1 2 0; 1 1 0 2];
b = [0.6; 0; 0; 0.4];

% call "cholesky()" to get R
R = cholesky(A);
fprintf("R = \n");
disp(R);
fprintf("R* = \n");
disp(R');
% verify A = R*R
fprintf("A = \n");
disp(R'*R);

% call "forward_backward_sub()" to get x
x = forward_backward_sub(R, b);

fprintf("x = \n");
disp(x);
