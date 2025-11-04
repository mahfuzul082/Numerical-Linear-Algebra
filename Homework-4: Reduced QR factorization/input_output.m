clear all; clc;
format longE;

% input matrix
A = [1 1+1e-8 1; 1 1 1+1e-8; 1 1 1];
fprintf("A = \n");
disp(A);

% call to classical gram-schmidt
[Q_GS, R_GS] = gram_schmidt(A);
fprintf("For Classical Gram-Schmidt:\n");
fprintf("Q = \n");
disp(Q_GS);
fprintf("R = \n");
disp(R_GS);

% call to modified gram-schmidt
[Q_MGS, R_MGS] = mod_gram_schmidt(A);
fprintf("For Modified Gram-Schmidt:\n");
fprintf("Q = \n");
disp(Q_MGS);
fprintf("R = \n");
disp(R_MGS);

% orthogonality check between 2nd and 3rd column
Ortho_23_GS = dot(Q_GS(:, 2), Q_GS(:, 3));
fprintf("Dot product value between 2nd and 3rd column for Classical Gram-Schmidt: %.16e\n", Ortho_23_GS);
if abs(Ortho_23_GS) < 1e-16
    fprintf("\n2nd and 3rd column are orthogonal to each other for Q in Classical Gram-Schmidt.\n");
else
    fprintf("\n2nd and 3rd column are not orthogonal to each other for Q in Classical Gram-Schmidt.\n");
end
Ortho_23_MGS = dot(Q_MGS(:, 2), Q_MGS(:, 3));
fprintf("\nDot product value between 2nd and 3rd column for Modified Gram-Schmidt: %.16e\n", Ortho_23_MGS);
if abs(Ortho_23_MGS) < 1e-16
    fprintf("\n2nd and 3rd column are orthogonal to each other for Q in Modified Gram-Schmidt.\n");
else
    fprintf("2nd and 3rd column are not orthogonal to each other for Q in Modified Gram-Schmidt.");
end